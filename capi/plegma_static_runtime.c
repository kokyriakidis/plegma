/* The Mojo runtime functions plegma's C API calls, so build/libplegma.a (pixi run build-static) needs no Mojo
 * runtime library and links into static executables.
 *
 * Mojo's runtime (libKGENCompilerRTShared.so, libAsyncRTMojoBindings.so) only ships as shared libraries.
 * plegma's C API uses 12 of its functions; this file implements them with libc and pthreads, following the
 * runtime's own sources (github.com/modular/modular: Mojo/lib/CompilerRT/{Memory,System,TracyBridge}.cpp,
 * MLRT/lib/Driver/DeviceContext/MojoBindings.cpp) and how max.algorithm.parallelize drives a CPU
 * DeviceContext (max/mojo/max/gpu/host/_legacy/device_context.mojo, enqueue_cpu_range).
 *
 * CPU tasks run on a lazy persistent pthread pool with:
 * - single slab allocation per enqueueHostFunctionRange (1 malloc/free per parallel region),
 * - targeted idle-worker wakeups (no thundering herd on small dispatches),
 * - own-context-first task popping in synchronize() (prevents cross-caller latency inversion while remaining
 *   deadlock-free under arbitrary nested parallelize() calls),
 * - pthread_atfork() safety and clean destructor shutdown for shared-library dlclose() / process exit.
 *
 * Signatures belong to the pinned Mojo nightly (pixi.toml). A runtime function plegma starts to call shows up as
 * an undefined symbol when test-static links; after moving the nightly, check these against its runtime sources and
 * run pixi run test-static. A Mojo program (a main) needs more of the runtime than this: it is for the library. */

#define _GNU_SOURCE
#include <pthread.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#ifdef __linux__
#include <sched.h>
#endif

/* ---- KGEN compiler runtime ---- */

/* Memory.cpp: alignment <= 0 means the preferred alignment; round up to a power of two >= sizeof(void*). */
void* KGEN_CompilerRT_AlignedAlloc(ssize_t alignment, ssize_t size)
{
    if (size < 0) {
        return NULL;
    }
    size_t a = alignment <= 0 ? 64 : (size_t)alignment;
    if (a < sizeof(void*)) {
        a = sizeof(void*);
    }
    /* Round non-power-of-two alignments up so posix_memalign never rejects with EINVAL. */
    if ((a & (a - 1)) != 0) {
        size_t pow2 = sizeof(void*);
        while (pow2 < a && pow2 <= (SIZE_MAX >> 1)) {
            pow2 <<= 1;
        }
        if (pow2 < a) {
            return NULL;
        }
        a = pow2;
    }
    void* p = NULL;
    if (posix_memalign(&p, a, size > 0 ? (size_t)size : 1) != 0) {
        return NULL;
    }
    return p;
}

void KGEN_CompilerRT_AlignedFree(void* p) { free(p); }

/* System.cpp: 0 means no stack trace (the runtime's default unless max-debug.stack-trace-on-error). */
int KGEN_CompilerRT_GetStackTrace(char** strings, unsigned depth)
{
    (void)strings;
    (void)depth;
    return 0;
}

int KGEN_CompilerRT_NumLogicalCores(void)
{
#ifdef __linux__
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        int count = CPU_COUNT(&set);
        if (count > 0) {
            return count;
        }
    }
#endif
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
}

/* TracyBridge.cpp, built without profiling. */
size_t KGEN_CompilerRT_TracyIsEnabled(void) { return 0; }
uint64_t KGEN_CompilerRT_TracyZoneBegin(const char* name, size_t length, uint32_t color)
{
    (void)name;
    (void)length;
    (void)color;
    return 0;
}
void KGEN_CompilerRT_TracyZoneEnd(uint64_t context) { (void)context; }

/* ---- AsyncRT CPU DeviceContext: persistent work-stealing pthread pool ----
 * parallelize: create("cpu"), deviceApi (must say "cpu"), enqueueHostFunctionRange with one coroutine per
 * work item (each resumed once runs its item to completion), synchronize, release. Errors are returned as
 * strings the caller frees with strfree; NULL is success. */

typedef void (*CoroutineFunction)(int8_t*);

typedef struct DeviceContext DeviceContext;
typedef struct TaskBatch TaskBatch;

typedef struct Task {
    CoroutineFunction resume;
    CoroutineFunction destroy;
    int8_t* handle;
    DeviceContext* owner;
    TaskBatch* batch;
    struct Task* next;
} Task;

struct TaskBatch {
    size_t pending_free;
    Task tasks[1];
};

struct DeviceContext {
    size_t remaining;
    pthread_cond_t done_cond;
};

typedef struct StringRef {
    const char* data;
    size_t length;
} StringRef;

#define PLEGMA_MAX_POOL_WORKERS 256

static const char kOomError[] = "plegma static runtime: out of memory";

static pthread_mutex_t g_init_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t g_pool_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t g_pool_cond = PTHREAD_COND_INITIALIZER;
static pthread_once_t g_atfork_once = PTHREAD_ONCE_INIT;
static Task* g_queue_head = NULL;
static Task* g_queue_tail = NULL;
static pthread_t g_workers[PLEGMA_MAX_POOL_WORKERS];
static int g_worker_count = 0;
static int g_idle_workers = 0;
static int g_pool_shutdown = 0;
static __thread int t_is_pool_worker = 0;

static char* error_string(const char* message)
{
    char* copy = strdup(message);
    return copy ? copy : (char*)kOomError;
}

static void atfork_prepare(void)
{
    pthread_mutex_lock(&g_init_mutex);
    pthread_mutex_lock(&g_pool_mutex);
}

static void atfork_parent(void)
{
    pthread_mutex_unlock(&g_pool_mutex);
    pthread_mutex_unlock(&g_init_mutex);
}

static void atfork_child(void)
{
    /* Reinitialize synchronization state in the child: parent threads do not exist across fork(). */
    pthread_mutex_init(&g_init_mutex, NULL);
    pthread_mutex_init(&g_pool_mutex, NULL);
    pthread_cond_init(&g_pool_cond, NULL);
    g_queue_head = NULL;
    g_queue_tail = NULL;
    g_worker_count = 0;
    g_idle_workers = 0;
    g_pool_shutdown = 0;
    t_is_pool_worker = 0;
}

static void register_atfork_once(void) { (void)pthread_atfork(atfork_prepare, atfork_parent, atfork_child); }

static Task* pop_any_locked(void)
{
    Task* task = g_queue_head;
    if (task) {
        g_queue_head = task->next;
        if (!g_queue_head) {
            g_queue_tail = NULL;
        }
        task->next = NULL;
    }
    return task;
}

static Task* pop_for_context_locked(DeviceContext* preferred)
{
    Task* prev = NULL;
    Task* cur = g_queue_head;
    while (cur) {
        if (cur->owner == preferred) {
            if (prev) {
                prev->next = cur->next;
            } else {
                g_queue_head = cur->next;
            }
            if (cur == g_queue_tail) {
                g_queue_tail = prev;
            }
            cur->next = NULL;
            return cur;
        }
        prev = cur;
        cur = cur->next;
    }
    return NULL;
}

static TaskBatch* complete_task_locked(DeviceContext* owner, TaskBatch* batch)
{
    owner->remaining--;
    if (owner->remaining == 0) {
        pthread_cond_broadcast(&owner->done_cond);
    }
    batch->pending_free--;
    return batch->pending_free == 0 ? batch : NULL;
}

static void* pool_worker(void* unused)
{
    (void)unused;
    t_is_pool_worker = 1;
    pthread_mutex_lock(&g_pool_mutex);
    for (;;) {
        while (!g_queue_head && !g_pool_shutdown) {
            g_idle_workers++;
            pthread_cond_wait(&g_pool_cond, &g_pool_mutex);
            g_idle_workers--;
        }
        if (g_pool_shutdown && !g_queue_head) {
            pthread_mutex_unlock(&g_pool_mutex);
            return NULL;
        }
        Task* task = pop_any_locked();
        pthread_mutex_unlock(&g_pool_mutex);

        DeviceContext* owner = task->owner;
        TaskBatch* batch = task->batch;
        task->resume(task->handle);
        task->destroy(task->handle);

        pthread_mutex_lock(&g_pool_mutex);
        TaskBatch* to_free = complete_task_locked(owner, batch);
        if (to_free) {
            pthread_mutex_unlock(&g_pool_mutex);
            free(to_free);
            pthread_mutex_lock(&g_pool_mutex);
        }
    }
}

static void ensure_pool_workers(size_t requested)
{
    pthread_once(&g_atfork_once, register_atfork_once);
    int max_workers = KGEN_CompilerRT_NumLogicalCores();
    if (max_workers < 1) {
        max_workers = 1;
    }
    if (max_workers > PLEGMA_MAX_POOL_WORKERS) {
        max_workers = PLEGMA_MAX_POOL_WORKERS;
    }
    int target = requested > (size_t)max_workers ? max_workers : (int)requested;
    if (target <= g_worker_count) {
        return;
    }
    /* Create threads under g_init_mutex (never holding g_pool_mutex across pthread_create). */
    pthread_mutex_lock(&g_init_mutex);
    if (!g_pool_shutdown && g_worker_count < target) {
        for (int i = g_worker_count; i < target; i++) {
            pthread_t tid;
            if (pthread_create(&tid, NULL, pool_worker, NULL) != 0) {
                break;
            }
            g_workers[g_worker_count++] = tid;
        }
    }
    pthread_mutex_unlock(&g_init_mutex);
}

__attribute__((destructor)) static void plegma_static_runtime_shutdown(void)
{
    pthread_mutex_lock(&g_init_mutex);
    pthread_mutex_lock(&g_pool_mutex);
    g_pool_shutdown = 1;
    pthread_cond_broadcast(&g_pool_cond);
    pthread_mutex_unlock(&g_pool_mutex);
    for (int i = 0; i < g_worker_count; i++) {
        (void)pthread_join(g_workers[i], NULL);
    }
    g_worker_count = 0;
    pthread_mutex_unlock(&g_init_mutex);
}

const char* AsyncRT_DeviceContext_create(DeviceContext** result, const char* kind, int id)
{
    (void)id;
    if (strcmp(kind, "cpu") != 0) {
        return error_string("plegma static runtime: only the CPU device is available");
    }
    DeviceContext* context = calloc(1, sizeof(DeviceContext));
    if (!context) {
        return error_string(kOomError);
    }
    if (pthread_cond_init(&context->done_cond, NULL) != 0) {
        free(context);
        return error_string("plegma static runtime: pthread_cond_init failed");
    }
    *result = context;
    return NULL;
}

void AsyncRT_DeviceContext_deviceApi(StringRef* result, DeviceContext* context)
{
    (void)context;
    result->data = "cpu";
    result->length = 3;
}

const char* AsyncRT_DeviceContext_enqueueHostFunctionRange(DeviceContext* context, CoroutineFunction resume,
                                                          CoroutineFunction destroy, int8_t** handles,
                                                          size_t count)
{
    if (count == 0) {
        return NULL;
    }
    if (count > (SIZE_MAX - sizeof(TaskBatch)) / sizeof(Task)) {
        for (size_t i = 0; i < count; i++) {
            destroy(handles[i]);
        }
        return error_string(kOomError);
    }
    TaskBatch* batch = malloc(sizeof(TaskBatch) + (count - 1) * sizeof(Task));
    if (!batch) {
        for (size_t i = 0; i < count; i++) {
            destroy(handles[i]);
        }
        return error_string(kOomError);
    }
    batch->pending_free = count;
    for (size_t i = 0; i < count; i++) {
        Task* task = &batch->tasks[i];
        task->resume = resume;
        task->destroy = destroy;
        task->handle = handles[i];
        task->owner = context;
        task->batch = batch;
        task->next = (i + 1 < count) ? &batch->tasks[i + 1] : NULL;
    }

    ensure_pool_workers(count);

    pthread_mutex_lock(&g_pool_mutex);
    if (g_queue_tail) {
        g_queue_tail->next = &batch->tasks[0];
    } else {
        g_queue_head = &batch->tasks[0];
    }
    g_queue_tail = &batch->tasks[count - 1];
    context->remaining += count;
    if (g_idle_workers > 0) {
        size_t wake = count < (size_t)g_idle_workers ? count : (size_t)g_idle_workers;
        if (wake == 1) {
            pthread_cond_signal(&g_pool_cond);
        } else if (wake == (size_t)g_idle_workers) {
            pthread_cond_broadcast(&g_pool_cond);
        } else {
            for (size_t w = 0; w < wake; w++) {
                pthread_cond_signal(&g_pool_cond);
            }
        }
    }
    pthread_mutex_unlock(&g_pool_mutex);
    return NULL;
}

const char* AsyncRT_DeviceContext_synchronize(DeviceContext* context)
{
    pthread_mutex_lock(&g_pool_mutex);
    while (context->remaining > 0) {
        /* Prefer our own context's tasks first so external callers avoid cross-caller latency inversion; pool
         * workers (nested spread) or fallback when no workers exist drain any task for deadlock freedom. */
        Task* task = pop_for_context_locked(context);
        if (!task && (t_is_pool_worker || g_worker_count == 0)) {
            task = pop_any_locked();
        }
        if (task) {
            pthread_mutex_unlock(&g_pool_mutex);
            DeviceContext* owner = task->owner;
            TaskBatch* batch = task->batch;
            task->resume(task->handle);
            task->destroy(task->handle);
            pthread_mutex_lock(&g_pool_mutex);
            TaskBatch* to_free = complete_task_locked(owner, batch);
            if (to_free) {
                pthread_mutex_unlock(&g_pool_mutex);
                free(to_free);
                pthread_mutex_lock(&g_pool_mutex);
            }
        } else {
            pthread_cond_wait(&context->done_cond, &g_pool_mutex);
        }
    }
    pthread_mutex_unlock(&g_pool_mutex);
    return NULL;
}

void AsyncRT_DeviceContext_release(DeviceContext* context)
{
    AsyncRT_DeviceContext_synchronize(context);
    pthread_cond_destroy(&context->done_cond);
    free(context);
}

void AsyncRT_DeviceContext_strfree(const char* p)
{
    if (p && p != kOomError) {
        free((void*)p);
    }
}

#if defined(__linux__) && defined(__GLIBC__) && (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
/* ---- glibc compatibility ----
 * The Mojo standard library calls __xstat and __lxstat, which glibc 2.33 and later keep only as shared-library
 * compatibility symbols, so a static link would not find them. _STAT_VER (1 on x86-64) uses struct stat. */
int __xstat(int version, const char* path, struct stat* buffer)
{
    (void)version;
    return stat(path, buffer);
}

int __lxstat(int version, const char* path, struct stat* buffer)
{
    (void)version;
    return lstat(path, buffer);
}
#endif
