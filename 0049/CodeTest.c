#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>
/*
#pragma warning(disable:6011)
#pragma warning(disable:6385)
#pragma warning(disable:6386)
*/

int main()
{
    int max, pnum = 0, pmax, i, j;
    char* pSieve;
    int* pPrime;
    /*
    printf("素数調査の最大数を入力してください (2以上の整数) \n");
    do
    {
        scanf_s("%d", &max);
    } while (max < 2);
    */
    max = 1000000000;
    printf("max = %d\n", max);
    max++;
    if (max < 10000)
    {
        pmax = max / 2;
    }
    else
    {
        pmax = max / log(max) * 1.2;
    }
    clock_t start = clock();
    pSieve = malloc(sizeof (char) * max / 2);
    if (!pSieve) exit(1);
    pPrime = malloc(sizeof(int) * pmax);
    if (!pPrime) exit(1);
    memset(pSieve, 1, max / 2);
    for (i =3; i * i < max; i += 2)
    {
        if (pSieve[i / 2])
        {
            for (j= i * i / 2; j< max / 2; j += i)
            {
                pSieve[j] = 0;
            }
        }
    }
    pPrime[0] = 2;
    pnum++;
    for (i = 1; i < max / 2; i++)
    {
        if (pSieve[i])
        {
            pPrime[pnum] = i*2+1;
            pnum++;
        }
    }
    clock_t t = clock() - start;
    
    printf("%ld.%ldsec\n", t / 1000, t % 1000);
    printf("素数が%d個見つかりました\n", pnum);
    free (pSieve);
    free (pPrime);
    return EXIT_SUCCESS;
}
