#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

void get_indices_3_ifo_forloop(int64_t tlen, int64_t t2_coinc_window, int64_t t3_coinc_window, int64_t* idx, int64_t num_combinations) {
    int64_t row_num = 0;
    for (int64_t i = 0; i < tlen; i++) {
        for (int64_t j = (i - t2_coinc_window > 0 ? i - t2_coinc_window : 0); j < tlen && j <= i + t2_coinc_window; j++) {
            for (int64_t k = (i - t3_coinc_window > 0 ? i - t3_coinc_window : 0); k < tlen && k <= i + t3_coinc_window; k++) {
                if (row_num < num_combinations) {
                    idx[row_num * 3] = i;
                    idx[row_num * 3 + 1] = j;
                    idx[row_num * 3 + 2] = k;
                    row_num++;
                } else {
                    return;
                }
            }
        }
    }
}

// int main() {
//     int64_t num_combinations = 35;
//     int64_t dtype = sizeof(int64_t);
//     int64_t tlen = 5;
//     int64_t t2_coinc_window = 1;
//     int64_t t3_coinc_window = 1;
//     int64_t *idx = malloc(num_combinations * 3 * dtype);
//     get_indices_3_ifo_forloop(tlen, t2_coinc_window, t3_coinc_window, idx, num_combinations);
//     for (int64_t i = 0; i < num_combinations; i++) {
//         printf("[%ld, %ld, %ld]\n", idx[i * 3], idx[i * 3 + 1], idx[i * 3 + 2]);
//     }
//     free(idx);
//     return 0;
// }
