## chemostat-hs

Implementing a chemostat model in Haskell, using [hmatrix](https://hackage.haskell.org/package/hmatrix).

### Dependencies

For building:
+ [stack](https://github.com/commercialhaskell/stack) ≥ 1.6.1

For plotting:

+ python3
+ matplotlib

### How to run

In the project root directory run `stack build` and then `stack exec chemostat`.

The plot file will be written to that directory as well.
