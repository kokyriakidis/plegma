/* pixi run test-static: parallelize through capi/plegma_static_runtime.c in a static executable. 64 items over 8
 * threads, 50 times, each time compared with one thread. Run with MALLOC_CHECK_=3 to catch a wrong free. */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

void plegma_test_parallel(int64_t* values, long items, long workers);

int main(void)
{
    enum { ITEMS = 64, REPEATS = 50 };
    int64_t serial[ITEMS], parallel[ITEMS];
    plegma_test_parallel(serial, ITEMS, 1);
    int differ = 0;
    for (int repeat = 0; repeat < REPEATS; repeat++) {
        memset(parallel, 0, sizeof parallel);
        plegma_test_parallel(parallel, ITEMS, 8);
        differ += memcmp(serial, parallel, sizeof serial) != 0;
    }
    printf("%s parallelize: %d of %d runs over 8 threads differ from one thread\n", differ ? "FAIL" : "PASS", differ,
           REPEATS);
    return differ != 0;
}
