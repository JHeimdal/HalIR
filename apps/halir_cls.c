#include "HalIR/halir.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

int main(void) {
    printf("HalIR version: %s\n", "1.0.0"); // Replace with actual version if available
    halir_simulation_setup *work = halir_simulation_setup_create();
    if (work == NULL) {
        fprintf(stderr, "Failed to create HalIR workspace\n");
        return 1;
    }
    halir_simulation_setup_free(work);
    return 0;
}