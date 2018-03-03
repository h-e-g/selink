# selink
Optimized window-based estimation of haplotype tests for positive selection

## Compile from source

If you've just clone the repository, you need to do once
```
$ aclocal
$ autoheader
$ autoconf
$ automake --add-missing
```
Then (or in any case), each time you want to recompile, do
```
$ ./configure
$ make
```
and optionnaly `make check`

If you need to regenerate the output of tests do
`make REGENERATE=true check`
