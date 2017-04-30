#include <stdio.h>
#include <stdlib.h>

int main() {

    FILE* fp;
    long double tok_int[3];

    fp = fopen("file.txt", "w+");

    fprintf(fp, "%.21Le %.21Le %.21Le \n", 0.3141592653589793116L, 3.141592653589793116L, 3.141592653589793116L);
    fprintf(fp, "%.21Le %.21Le %.21Le \n", 0.3141592653589793116L, 3.141592653589793116L, 3.141592653589793116L);


    fclose(fp);

    fp = fopen("file.txt", "r");

    while(!feof(fp)) {
        fscanf(fp,"%Lf %Lf %Lf \n", &tok_int[0], &tok_int[1], &tok_int[2]);
        printf("%.18Le \n", tok_int[2]);
    }

    fclose(fp);

    return 0;
}