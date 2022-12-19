#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hitran.h"

void readParameters(char *filename) {
    HitranHead head;
    HitranLine *lines;
    FILE *prm_file = fopen(filename, "rb");
    fread( (void*)&head, sizeof(HitranHead), 1, prm_file);
    head.molecs = (int*)malloc(sizeof(int)*head.nisotp);
    fread( (void*)head.molecs, sizeof(int), head.nisotp, prm_file);

    lines = (HitranLine*)malloc(sizeof(HitranLine)*head.ndatapnts);
    for (int i=0;i<head.ndatapnts;i++)
        fread( &lines[i], sizeof(HitranLine), 1, prm_file);

    printf("%s\n", head.molec);
    printf("%s\n", head.isoname);
    printf("%d\n", head.nisotp);
    printf("%f\n", head.roi_low);
    printf("%f\n", head.roi_high);
    printf("%d\n", head.ndatapnts);
    printf("%d\n", head.molecs[0]);
    for (int i=0;i<head.ndatapnts;i++) {
        printf("%d %d %f %e\n", lines[i].molec_num, lines[i].isotp_num, lines[i].trans_mu, lines[i].line_I);
    }
}

void calcSpectra(void) {
}


int main()
{
    readParameters("CO.hpar");
    return 0;
}
