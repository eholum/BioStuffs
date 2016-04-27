/*
 * Counts lines/character occurences in a file. Usefull for fasta/fastq files.
 * 
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#define ASCII_SIZE 256 /* Extended ascii */

typedef unsigned long count; /* Counter */

/* Current file counters */
static count ccount[ASCII_SIZE];
count lcount;

/* Total files counter */
static count t_ccount[ASCII_SIZE];
count t_lcount;

void report (char *file, count *ccount, count lcount, char *s) 
{
    int i;
    for (i = 0; *(s + i) != '\0'; i++) {
        printf("%6lu ", ccount[*(s + i)]);
    }
    printf("%6lu %s\n ", lcount, file);
}

void reset (count *ccount, char *s) 
{
    int i;
    memset(ccount, -1L, ASCII_SIZE);
    for (i = 0; *(s + i) != '\0'; i++) {
       ccount[*(s + i)] = 0L;    
    }
}

void counter(char *file, char *s) 
{
    unsigned int c;
    
    FILE *fp = fopen(file, "r");
    
    if (!fp) {
        printf ("cannot open file '%s'\n", file);
        return;
    }
    
    reset(ccount, s);
    lcount = 0L;
    while ((c = getc(fp)) != EOF ) {
        // Increment line count first
        if (c == '\n') {
            lcount++; 
            t_lcount++;
        }
        // If we don't care about the char, skip it
        if (ccount[c] == -1L)
            continue;
        
        (ccount[c])++;
        (t_ccount[c])++;
    }
    
    report (file, ccount, lcount, s);
}
    
int main(int argc, char **argv) 
{
    int ch, i;
    
    if (argc < 3) {
        printf ("usage: ccount CHARS FILE [FILE...]\n");
        return -1;
    }
    
    for (i = 0; argv[1][i] != '\0'; i++) {
        printf("%6c ", argv[1][i]);
    }
    printf("lines file \n");
    
    reset (t_ccount, argv[1]);
    for (i = 2; i < argc; i++)
        counter (argv[i], argv[1]);
    
    if (argc > 3)
        report ("total", t_ccount, t_lcount, argv[1]);
    
    return 0;
}