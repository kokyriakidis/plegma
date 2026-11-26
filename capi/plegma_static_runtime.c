/* The Mojo runtime functions plegma's C API calls, so build/libplegma.a (pixi run build-static) needs no Mojo
 * runtime library and links into static executables.
 *
 * Mojo's runtime (libKGENCompilerRTShared.so, libAsyncRTMojoBindings.so) only ships as shared libraries.
 * plegma's C API uses 12 of its functions; this file implements them with libc and pthreads, following the
 * runtime's own sources (github.com/modular/modular: Mojo/lib/CompilerRT/{Memory,System,TracyBridge}.cpp,
 * MLRT/lib/Driver/DeviceContext/MojoBindings.cpp) and how max.algorithm.parallelize drives a CPU
 * DeviceContext (max/mojo/max/gpu/host/_legacy/device_context.mojo, enqueue_cpu_range).
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

/* ---- AsyncRT CPU DeviceContext, as parallelize uses it ----
 * parallelize: create("cpu"), deviceApi (must say "cpu"), enqueueHostFunctionRange with one coroutine per
 * work item (each resumed once runs its item to completion), synchronize, release. Errors are returned as
 * strings the caller frees with strfree; NULL is success. */

typedef void (*CoroutineFunction)(int8_t*);

typedef struct Task {
    CoroutineFunction resume;
    CoroutineFunction destroy;
    int8_t* handle;
} Task;

typedef struct DeviceContext {
    pthread_mutex_t mutex;
    pthread_t* threads;
    int* joinable;
    size_t count;
    size_t capacity;
} DeviceContext;

typedef struct StringRef {
    const char* data;
    size_t length;
} StringRef;

static char* error_string(const char* message) { return strdup(message); }

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
    pthread_mutex_init(&context->mutex, NULL);
    *result = context;
    return NULL;
}

void AsyncRT_DeviceContext_deviceApi(StringRef* result, DeviceContext* context)
{
    (void)context;
    result->data = "cpu";
    result->length = 3;
}

/* Runs one coroutine to completion and destroys it; frees its Task. */
static void* run_task(void* argument)
{
    Task* task = argument;
    task->resume(task->handle);
    task->destroy(task->handle);
    free(task);
    return NULL;
}

/* Each coroutine on a thread of its own: parallelize makes one per worker, and karan asks for no more
 * workers than the machine runs at once. */
const char* AsyncRT_DeviceContext_enqueueHostFunctionRange(DeviceContext* context, CoroutineFunction resume,
                                                          CoroutineFunction destroy, int8_t** handles,
                                                          size_t count)
{
    pthread_mutex_lock(&context->mutex);
    if (context->count + count > context->capacity) {
        size_t capacity = (context->count + count) * 2;
        pthread_t* threads = realloc(context->threads, capacity * sizeof(pthread_t));
        if (threads) {
            context->threads = threads;
        }
        int* joinable = realloc(context->joinable, capacity * sizeof(int));
        if (joinable) {
            context->joinable = joinable;
        }
        if (!threads || !joinable) {
            pthread_mutex_unlock(&context->mutex);
            return error_string("plegma static runtime: out of memory");
        }
        context->capacity = capacity;
    }
    for (size_t i = 0; i < count; i++) {
        Task* task = malloc(sizeof(Task));
        if (!task) {
            for (size_t j = i; j < count; j++) {
                destroy(handles[j]);
            }
            pthread_mutex_unlock(&context->mutex);
            return error_string("plegma static runtime: out of memory");
        }
        task->resume = resume;
        task->destroy = destroy;
        task->handle = handles[i];
        context->joinable[context->count] =
            pthread_create(&context->threads[context->count], NULL, run_task, task) == 0;
        if (!context->joinable[context->count]) {
            /* No thread: run it here. */
            run_task(task);
        }
        context->count++;
    }
    pthread_mutex_unlock(&context->mutex);
    return NULL;
}

const char* AsyncRT_DeviceContext_synchronize(DeviceContext* context)
{
    pthread_mutex_lock(&context->mutex);
    for (size_t i = 0; i < context->count; i++) {
        if (context->joinable[i]) {
            pthread_join(context->threads[i], NULL);
        }
    }
    context->count = 0;
    pthread_mutex_unlock(&context->mutex);
    return NULL;
}

void AsyncRT_DeviceContext_release(DeviceContext* context)
{
    AsyncRT_DeviceContext_synchronize(context);
    pthread_mutex_destroy(&context->mutex);
    free(context->joinable);
    free(context->threads);
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
