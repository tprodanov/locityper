// NOTE: Compile with `gcc rename_fq.c -o bin/rename_fq -O3`

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main (int argc, char *argv[]) {
    char prefix[128];
    assert(argc == 3);
    strcpy(prefix, basename(argv[1]));
    size_t len = strlen(prefix);
    assert(len > 5);
    prefix[len - 5] = 0;

    FILE * f1 = fopen(argv[1], "r");
    FILE * f2 = fopen(argv[2], "r");
    assert(f1 != NULL && f2 != NULL);

    size_t size = 0;
    size_t read;
    char* line;
    size_t num = 1;
    while (getline(&line, &size, f1) > 0) {
        printf("@%s:%X/1\n", prefix, num);
        assert(getline(&line, &size, f1) != -1);
        printf(line);
        assert(getline(&line, &size, f1) != -1);
        printf("+\n");
        assert(getline(&line, &size, f1) != -1);
        printf(line);

        assert(getline(&line, &size, f2) != -1);
        printf("@%s:%X/2\n", prefix, num);
        assert(getline(&line, &size, f2) != -1);
        printf(line);
        assert(getline(&line, &size, f2) != -1);
        printf("+\n");
        assert(getline(&line, &size, f2) != -1);
        printf(line);
        ++num;
    }

    fclose(f1);
    fclose(f2);
    if (line)
        free(line);
    exit(EXIT_SUCCESS);
}
