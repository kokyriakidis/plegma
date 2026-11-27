/* The Mojo runtime functions plegma's C API calls, so build/libplegma.a (pixi run build-static) needs no Mojo
 * runtime library and links into static executables.
 *
 * Mojo's runtime (libKGENCompilerRTShared.so, libAsyncRTMojoBindings.so) only ships as shared libraries.
 * plegma's C API uses 12 of its functions; this file implements them with libc and pthreads, following the
 * runtime's own sources (github.com/modular/modular: Mojo/lib/CompilerRT/{Memory,System,TracyBridge}.cpp,
 * MLRT/lib/Driver/DeviceContext/MojoBindings.cpp) and how max.algorithm.parallelize drives a CPU
 * DeviceContext (max/mojo/max/gpu/host/_legacy/device_context.mojo, enqueue_cpu_range).
 *
 * CPU tasks run on a lazy persistent pthread pool with caller-helps-drain synchronization, so nested `spread`
 * calls never deadlock and repeated parallel regions avoid per-call pthread creation overhead.
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

/* ---- KGEN compiler runtime ---- */

/* Memory.cpp: alignment <= 0 means the preferred alignment. */
void* KGEN_CompilerRT_AlignedAlloc(ssize_t alignment, ssize_t size)
{
    size_t a = alignment <= 0 ? 64 : (size_t)alignment;
    if (a < sizeof(void*)) {
        a = sizeof(void*);
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
 * strings the caller frees with strfree; NULL is success.
 *
 * A single mutex guards the shared work queue and every DeviceContext's remaining counter. Waiting threads in
 * synchronize() pop and execute queued tasks while waiting, guaranteeing deadlock-freedom under arbitrary nested
 * parallelize() calls and concurrent external caller threads. */

typedef void (*CoroutineFunction)(int8_t*);

typedef struct DeviceContext DeviceContext;

typedef struct Task {
    CoroutineFunction resume;
    CoroutineFunction destroy;
    int8_t* handle;
    DeviceContext* owner;
    struct Task* next;
} Task;

struct DeviceContext {
    size_t remaining;
    pthread_cond_t done_cond;
};

typedef struct StringRef {
    const char* data;
    size_t length;
} StringRef;

static pthread_mutex_t g_pool_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t g_pool_cond = PTHREAD_COND_INITIALIZER;
static Task* g_queue_head = NULL;
static Task* g_queue_tail = NULL;
static int g_pool_started = 0;

static char* error_string(const char* message) { return strdup(message); }

static Task* pop_task_locked(void)
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

static void complete_task_locked(DeviceContext* owner)
{
    owner->remaining--;
    if (owner->remaining == 0) {
        pthread_cond_broadcast(&owner->done_cond);
    }
}

static void* pool_worker(void* unused)
{
    (void)unused;
    pthread_mutex_lock(&g_pool_mutex);
    for (;;) {
        while (!g_queue_head) {
            pthread_cond_wait(&g_pool_cond, &g_pool_mutex);
        }
        Task* task = pop_task_locked();
        pthread_mutex_unlock(&g_pool_mutex);

        DeviceContext* owner = task->owner;
        task->resume(task->handle);
        task->destroy(task->handle);
        free(task);

        pthread_mutex_lock(&g_pool_mutex);
        complete_task_locked(owner);
    }
    return NULL;
}

static void ensure_pool_started_locked(void)
{
    if (g_pool_started) {
        return;
    }
    g_pool_started = 1;
    int workers = KGEN_CompilerRT_NumLogicalCores();
    if (workers < 1) {
        workers = 1;
    }
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    for (int i = 0; i < workers; i++) {
        pthread_t tid;
        (void)pthread_create(&tid, &attr, pool_worker, NULL);
    }
    pthread_attr_destroy(&attr);
}

const char* AsyncRT_DeviceContext_create(DeviceContext** result, const char* kind, int id)
{
    (void)id;
    if (strcmp(kind, "cpu") != 0) {
        return error_string("plegma static runtime: only the CPU device is available");
    }
    DeviceContext* context = calloc(1, sizeof(DeviceContext));
    if (!context) {
        return error_string("plegma static runtime: out of memory");
    }
    pthread_cond_init(&context->done_cond, NULL);
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
    Task* head = NULL;
    Task* tail = NULL;
    for (size_t i = 0; i < count; i++) {
        Task* task = malloc(sizeof(Task));
        if (!task) {
            while (head) {
                Task* next = head->next;
                destroy(head->handle);
                free(head);
                head = next;
            }
            for (size_t j = i; j < count; j++) {
                destroy(handles[j]);
            }
            return error_string("plegma static runtime: out of memory");
        }
        task->resume = resume;
        task->destroy = destroy;
        task->handle = handles[i];
        task->owner = context;
        task->next = NULL;
        if (tail) {
            tail->next = task;
        } else {
            head = task;
        }
        tail = task;
    }
    pthread_mutex_lock(&g_pool_mutex);
    ensure_pool_started_locked();
    if (g_queue_tail) {
        g_queue_tail->next = head;
    } else {
        g_queue_head = head;
    }
    g_queue_tail = tail;
    context->remaining += count;
    if (count == 1) {
        pthread_cond_signal(&g_pool_cond);
    } else {
        pthread_cond_broadcast(&g_pool_cond);
    }
    pthread_mutex_unlock(&g_pool_mutex);
    return NULL;
}

const char* AsyncRT_DeviceContext_synchronize(DeviceContext* context)
{
    pthread_mutex_lock(&g_pool_mutex);
    while (context->remaining > 0) {
        Task* task = pop_task_locked();
        if (task) {
            pthread_mutex_unlock(&g_pool_mutex);
            DeviceContext* owner = task->owner;
            task->resume(task->handle);
            task->destroy(task->handle);
            free(task);
            pthread_mutex_lock(&g_pool_mutex);
            complete_task_locked(owner);
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

void AsyncRT_DeviceContext_strfree(const char* p) { free((void*)p); }

#ifdef __linux__
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
