#include <stdio.h>
int main() {
    const int num = 35;
    const unsigned long arr[] = {
        num, 1, 4, 3, 5, 14,
        2, 8, 9, 10, 4,
        22, 18, 7, 5, 34,
        19, 20, 17, 10, 7,
        8, 25, 27, 29, 13,
        15, 3, 2, 17, 15,
        24, 25, 28, 27, 20
    };
    FILE* file;
    file = fopen("binary_input", "wb");
    fwrite(arr, sizeof(arr), 1, file);
}