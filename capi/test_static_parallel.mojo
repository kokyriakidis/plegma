# For pixi run test-static: Mojo's parallelize (through karan's spread), called from C, so the CPU device functions
# in capi/plegma_static_runtime.c run. Each item allocates and computes; the C side checks threads give what one
# thread gives.
from karan.common import spread

comptime Out = Pointer[Int64, MutUntrackedOrigin]


@export
def plegma_test_parallel(values: Out, items: Int, workers: Int) abi("C"):
    def work(item: Int) {imm values}:
        var scratch = List[Int64](length=1000 + item, fill=0)
        var total: Int64 = 0
        for k in range(200000):
            scratch[k % len(scratch)] += Int64((item * 31 + k) % 7)
        for k in range(len(scratch)):
            total += scratch[k]
        values[unsafe_offset=item] = total

    spread(work, items, workers)
