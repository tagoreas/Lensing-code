#include <stdio.h>
#include <stdlib.h>

extern "C"
void ps_loaddata(char *a, size_t b, void *c, void *d)
{
    printf("This is a dummy version of pixsrc. Recompile for real functionality.\n");
    exit(1);
}
extern "C"
void ps_getpixsrcinfo(char *a, size_t *b, size_t *c, size_t *d)
{
    printf("This is a dummy version of pixsrc. Recompile for real functionality.\n");
    exit(1);
}
extern "C"
void ps_launch(void *a, void *b, double *c)
{
    printf("This is a dummy version of pixsrc. Recompile for real functionality.\n");
    exit(1);
}
extern "C"
void ps_freepixsrcmem(void *a, void *b)
{
    printf("This is a dummy version of pixsrc. Recompile for real functionality.\n");
    exit(1);
}
