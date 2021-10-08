#include <iostream>
#include <omp.h>
#include <unistd.h>

void first_task() {
#pragma omp parallel num_threads(8)
    {
        int count = omp_get_num_threads();
        int current_thread_num = omp_get_thread_num();
        printf("Hello world! current - %d, all - %d \n", current_thread_num, count);
    }
}
void second_task() {
    const int min = 2;
    int threads_count = 3;
#pragma omp parallel num_threads(threads_count) if(threads_count > min)
    {
        printf("first, %d out of %d \n", omp_get_thread_num(), omp_get_num_threads());
    }
    threads_count = 2;
#pragma omp parallel num_threads(threads_count) if (threads_count > min)
    {
        printf("second, %d out of %d \n", omp_get_thread_num(), omp_get_num_threads());
    }
}
void third_task() {
    int a = 0;
    int b = 0;
    printf("before first block a - %d, b - %d \n", a, b);
#pragma omp parallel num_threads(2) private(a) firstprivate(b)
    {
        int current_thread_num = omp_get_thread_num();
        a = 0;
        a += current_thread_num;
        b += current_thread_num;
        printf("current_thread_num - %d, a - %d, b - %d \n", current_thread_num, a, b);
    }
    printf("after first block a - %d, b - %d \n", a, b);
    printf("before second block a - %d, b - %d \n", a, b);
#pragma omp parallel num_threads(4) shared(a) private(b)
    {
        b = 4;
        int current_thread_num = omp_get_thread_num();
        a -= current_thread_num;
        b -= current_thread_num;
        printf("current_thread_num - %d, a - %d, b - %d \n", current_thread_num, a, b);
    }
    printf("after second block a - %d, b - %d \n", a, b);
}
void fouth_task()
{
    int a[10];
    int b[10];
    for (int i = 0; i < 10; i++)
    {
        a[i] = rand();
        b[i] = rand();
    }
#pragma omp parallel num_threads(2)
    {
#pragma omp master
        {
            int min = a[0];
            for (int i = 1; i < 10; i++)
            {
                if (a[i] < min) min = a[i];
            }
            printf("min a = %d, main thread = %d \n", min, omp_get_thread_num());
        }
#pragma omp single
        {
            int max = b[0];
            for (int i = 1; i < 10; i++)
            {
                if (b[i] > max) max = a[i];
            }
            printf("max b = %d, slave thread = %d \n", max, omp_get_thread_num());
        }
    }
}
void fifth_task()
{
    int d[6][8];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 8; j++) {
            d[i][j] = rand();
        }
    }
#pragma omp parallel sections num_threads(3)
    {
#pragma omp section
        {
            int sum = 0;
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 8; j++) {
                    sum += d[i][j];
                }
            }
            printf("first section: %d, thread_num = %d \n", sum / (6 * 8), omp_get_thread_num());
        }
#pragma omp section
        {
            int min = d[0][0];
            int max = d[0][0];
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 8; j++) {
                    if (d[i][j] > max) max = d[i][j];
                    if (d[i][j] < min) min = d[i][j];
                }
            }
            printf("min = %d, max = %d, thread_num = %d \n", min, max, omp_get_thread_num());
        }
#pragma omp section
        {
            int count = 0;
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 8; j++) {
                    if (d[i][j] % 3 == 0) count++;
                }
            }
            printf("devide by 3 count = %d, thread - %d \n", count, omp_get_thread_num());
        }
    }
}
void sixth_task()
{
    int a[100];
    for (int i = 0; i < 100; i++)
    {
        a[i] = rand();
    }
    int fsum = 0;
#pragma omp parallel for
    for(int i = 0; i < 100; i++)
    {
        fsum += a[i];
    }
    printf("for sum : %d \n", fsum);
    fsum = 0;
#pragma omp parallel for reduction(+:fsum)
    for (int i = 0; i < 100; i++) {
        {
            fsum += a[i];
        }
    }
    printf("for reduction sum : %d \n", fsum);
    fsum = 0;
    for (int i = 0; i < 100; i++)
    {
        fsum += a[i];
    }
    printf("for sum : %d \n", fsum);
}

void seventh_task()
{
    int a[12];
    int b[12];
    int c[12];
#pragma omp parallel for schedule(static, 4) num_threads(3)
    for (int i = 0; i < 12; i++)
    {
        a[i] = rand();
        b[i] = rand();
        printf("a[%d] = %d, b[%d] = %d, thread = %d \n", i, a[i], i, b[i], omp_get_thread_num());
    }
#pragma omp parallel for schedule(dynamic, 3) num_threads(4)
    for (int i = 0; i < 12; i++)
    {
        c[i] = a[i] + b[i];
        printf("c[%d]  = %d, a[%d] = %d, b[%d] = %d, thread = %d \n",i, c[i], i, a[i], i, b[i], omp_get_thread_num());
    }
}

void eigth_task_static(long a[], double b[]) {
    double start_time = omp_get_wtime();
#pragma omp parallel for schedule(static, 2000)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = ( a[i - 1] + a[i] + a[i + 1]) / 3.0;
    }
    double end = omp_get_wtime();
    printf("static time = %.16g \n", end - start_time);
}

void eigth_task_dynamic(long a[], double b[]) {
    double start_time = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 500)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
    }
    double end = omp_get_wtime();
    printf("dynamic time = %.16g \n", end - start_time);
}

void eigth_task_guided(long a[], double b[]) {
    double start_time = omp_get_wtime();
#pragma omp parallel for schedule(guided, 500)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
    }
    double end = omp_get_wtime();
    printf("guided time = %.16g \n", end - start_time);
}

void eigth_task_runtime(long a[], double b[]) {
    double start_time = omp_get_wtime();
#pragma omp parallel for schedule(runtime)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = (a[i - 1] + a[i] + a[i + 1]) / 3.0;
    }
    double end = omp_get_wtime();
    printf("runtime time = %.16g \n", end - start_time);
}

void eigth_task()
{
    long a[16000];
    double b[16000];
#pragma omp parallel for
    for (int i = 0; i < 16000; i++) {
        a[i] = i;
    }
    eigth_task_runtime(a, b);
    eigth_task_static(a, b);
    eigth_task_dynamic(a, b);
    eigth_task_guided(a, b);
}
void ninth_task()
{
    const int size = 2000;
    int* vector = new int[size];
    int** matrix = new int*[size];
    for (int i = 0; i < size; i++)
    {
        matrix[i] = new int[size];
        vector[i] = rand();
        for (int j = 0; j < size; j++)
        {
            matrix[i][j] = rand();
        }
    }
    double start_time = omp_get_wtime();
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            vector[i] = matrix[i][j] * vector[j];
        }
    }
    double end_time = omp_get_wtime();
    printf("seq time =  %.16g \n", end_time - start_time);
    start_time = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 63)
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            vector[i] = matrix[i][j] * vector[j];
        }
    }
    end_time = omp_get_wtime();
    printf("par time =  %.16g \n", end_time - start_time);
    delete[] vector;
    for (int i = 0; i < size; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;
}
void tenth_task()
{
    int d[6][8];
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            d[i][j] = rand();
        }
    }
    int max = d[0][0];
    int min = d[0][0];
#pragma omp parallel for
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            if (d[i][j] >= min && d[i][j] <= max) {
                continue;
            }
#pragma omp critical
            {
                if (d[i][j] < min) min = d[i][j];
                if (d[i][j] > max) max = d[i][j];
            }
        }
    }
    printf("max = %d, min = %d \n", min, max);
}
void eleventh_task()
{
    int a[30];
    int count = 0;
    for (int i = 0; i < 30; i++)
    {
        a[i] = rand();
    }
#pragma omp parallel for num_threads(8)
    for (int i = 0; i < 30; i++)
    {
        if (a[i] % 9 == 0)
        {
#pragma omp atomic
            count++;
        }
    }
    printf("divided by 9 count = %d", count);
}
void twelveth_task()
{
    int a[30];
    for (int i = 0; i < 30; i++)
    {
        a[i] = rand();
    }
    int max = INT16_MIN;
#pragma omp parallel for
    for (int i = 0; i < 30; i++)
    {
        if (a[i] % 7 == 0)
        {
#pragma omp critical
            {
                if (max < a[i]) max = a[i];
            }
        }
    }
    printf("max divided by 7 = %d", max);
}
void thirteenth_task_1()
{
    int current_thread = 7;
#pragma omp parallel num_threads(8)
    {
        bool isDone = false;
        while (!isDone)
        {
#pragma omp critical
            {
                if (omp_get_thread_num() == current_thread)
                {
                    current_thread--;
                }
            }
            printf("thread = %d \n", omp_get_thread_num());
            isDone = true;
        }
    }
}

void thirteenth_task_2()
{
#pragma omp parallel num_threads(8)
    {
        sleep(1000 / (omp_get_thread_num() + 1));
        printf("Thread %d \n", omp_get_thread_num());
    }
}

void thirteenth_task_3()
{
#pragma omp parallel
    {
        int nthreads = omp_get_num_threads();
        for (int i = nthreads - 1; i >= 0; i--)
        {
#pragma omp barrier
            {
                if (i == omp_get_thread_num())
                {
                    printf("thread = %d \n", omp_get_thread_num());
                }
            }
        }
    }
}

void thirteenth_task_4()
{
#pragma omp parallel for schedule(static, 1)
    for (int i = omp_get_num_threads() - 1; i >=0; i--)
    {
        sleep(100*i);
        printf("thread = %d \n", omp_get_thread_num());
    }
}

int thirteenth_task_5() {
    omp_lock_t my_lock;
    omp_init_lock(&my_lock);
    int current_thread = 7;
#pragma omp parallel
    {
        while (current_thread > 0)
        {
            omp_set_lock(&my_lock);
            if (current_thread == omp_get_thread_num()) {
                printf("thread = %d \n", omp_get_thread_num());
                current_thread--;
            }
            omp_unset_lock(&my_lock);
        }
    }
    omp_destroy_lock(&my_lock);
}

int main() {
    thirteenth_task_3();
    return EXIT_SUCCESS;
}
