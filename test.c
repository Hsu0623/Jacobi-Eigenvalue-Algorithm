#include <stdio.h>
#include <stdlib.h>

int main(){
    FILE *fp;
    int i;
    char filename[20] = "abc00";
    for(i = 0; i < 3; i++){
        fp = fopen(filename, "w");
        fprintf(fp, "1231");
        fclose(fp);
        if(filename[4]=='9'){
            filename[4] = 0;
            filename[3] += 1;
        }
        else filename[4]+=1;
    }
    return 0;
}