#include <stdio.h>
#include <stdlib.h>
int main() {
    unsigned long times;
    FILE* file;
    file = fopen("binary_input", "rb");
    fread(&times, sizeof(unsigned long), 1, file);
    fseek(file, sizeof(unsigned long), SEEK_SET);
    unsigned long* arr = (unsigned long*) malloc(sizeof(unsigned long) * times);
    fread(arr, sizeof(unsigned long) * times, 1, file);
    for (int i = 0; i < times; ++i) {
        printf("%lu\n", arr[i]);
    }
}