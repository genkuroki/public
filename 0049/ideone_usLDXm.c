#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>

int main()
{
    int max, pnum = 0, pmax, i, j;
    char* pSieve;
    int* pPrime;
    printf("素数調査の最大数を入力してください（2以上の整数）\n");
    do
    {
        scanf("%d", &max);
    } while (max < 2);
    clock_t start = clock();
    max++;
    if (max < 10000)
    {
        pmax = max / 2+1;
    }
    else
    {
        pmax = max / log(max) * 1.2;
    }
    pSieve = malloc(sizeof(char) * max / 2 + 16);
    if (!pSieve) exit(1);
    pPrime = malloc(sizeof(int) * pmax);
    if (!pPrime) exit(1);
    int* ipSieve = pSieve;
    for (j = 0; j <= max / 8 ; j += 3)
    {
        ipSieve[j + 0] = 0x01010001;
        ipSieve[j + 1] = 0x00010100;
        ipSieve[j + 2] = 0x01000101;
    }
    pSieve[1] = 1;
    clock_t lap1 = clock() - start;
    for (i = 5; i * i < max; i+=2)
    {
        if (pSieve[i / 2])
        {
            j = i * i / 2;
            pSieve[j] = 0;
            j += i * (i % 3 - 1);
            for (; j < max / 2 - i * 2; j += i)
            {
                pSieve[j] = 0;
                j += i * 2;
                pSieve[j] = 0;
            }
            if (j < max / 2)
            {
                pSieve[j] = 0;
            }
        }
    }
    clock_t lap2 = clock() - start;
    pPrime[0] = 2;
    pnum++;
    for (i = 1; i < max / 2; i++)
    {
        pPrime[pnum] = i*2+1;
        pnum+= pSieve[i];
    }
    clock_t lap3 = clock() - start;
    
    printf("%ld.%03ldsec\n", lap1 / 1000, lap1 % 1000);
    printf("%ld.%03ldsec\n", lap2 / 1000, lap2 % 1000);
    printf("%ld.%03ldsec\n", lap3 / 1000, lap3 % 1000);
    printf("素数が%d個見つかりました\n", pnum);
    free(pSieve);
    free(pPrime);
    return EXIT_SUCCESS;
}
