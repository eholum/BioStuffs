/*
 * Counts lines/character occurences in a file. Usefull for fasta/fastq files.
 * 
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

typedef unsigned long count; /* Counter */

/* Current file counters */
static count ccount[256];
count lcount;

/* Total files counter */
static count t_ccount[256];
count t_lcount;

void report (char *file, count *ccount, count lcount, char *s) 
{
    int i;
    
    for (i = 0; s[i] != '\0'; i++) {
        printf("%6c ", *(s + i));
    }
    printf("lines file \n");
    for (i = 0; *(s + i) != '\0'; i++) {
        printf("%6lu ", ccount[*(s + i)]);
    }
    printf("%6lu %s\n ", lcount, file);
}

void counter(char *file, char *s) 
{
    char c;
    
    FILE *fp = fopen(file, "r");
    
    if (!fp) {
        printf ("cannot open file '%s'", file);
        return;
    }
    
    memset(ccount, 0L, 256);
    while ((c = getc(fp)) != EOF) {
        if (c == '\n') {
            lcount++; t_lcount++;
        }
        (ccount[c])++; (t_ccount[c])++;
    }
    
    report (file, ccount, lcount, s);
}
    
int main(int argc, char **argv) 
{
    int ch, i;
    
    if (argc < 3) {
        printf ("usage: ccount CHARS FILE [FILE...]");
        return -1;
    }
    
    for (i = 2; i < argc; i++)
        counter (argv[i], argv[1]);
    
    if (argc > 3)
        report ("total", t_ccount, t_lcount, argv[1]);
    
    return 0;
}