#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;


void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}


int partition(int arr[], int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] <= pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}


int partition_r(int arr[], int low, int high) {
    srand(time(NULL));
    int random = low + rand() % (high - low);
    swap(&arr[random], &arr[high]);
    return partition(arr, low, high);
}


int randomizedSelect(int arr[], int low, int high, int k) {
    if (low == high)
        return arr[low];

    int pi = partition_r(arr, low, high);
    int length = pi - low + 1;

    if (length == k)
        return arr[pi];
    else if (k < length)
        return randomizedSelect(arr, low, pi - 1, k);
    else
        return randomizedSelect(arr, pi + 1, high, k - length);
}

int main() {
    int n, k;
    cout << "Enter the size of the array: ";
    cin >> n;
    int arr[n];
    cout << "Enter the elements of the array: ";
    for (int i = 0; i < n; i++) {
        cin >> arr[i];
    }
    cout << "Enter the value of k: ";
    cin >> k;

    int result = randomizedSelect(arr, 0, n - 1, k);
    cout << "The " << k << "-th smallest element is: " << result << endl;

    return 0;
}

