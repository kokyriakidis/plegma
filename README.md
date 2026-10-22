# plegma

Fast partial order alignment (POA) of DNA in Mojo: multiple sequence alignment and consensus from a graph that
each sequence is woven into. The goal is to be the fastest POA aligner by a wide margin; rivals are
[abPOA](https://github.com/yangao07/abPOA), [SPOA](https://github.com/rvaser/spoa),
[POASTA](https://github.com/broadinstitute/poasta) and [Theseus](https://github.com/Marco-Sola-Lab/theseus-lib).

*πλέγμα* (plegma) is Greek for a woven mesh or lattice, and is what a POA graph is: sequences threaded into one
another where they agree.

Status: early. The package holds the DNA alphabet; the graph, the alignment kernels and the command-line tool come
next.

## Building

plegma builds with [pixi](https://pixi.sh), which fetches the pinned Mojo nightly:

```sh
pixi run test     # build and run the tests
pixi run build    # precompile the package to build/plegma.mojoc
pixi run format   # mojo format, 120 columns
```

## License

MIT
