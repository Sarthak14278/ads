#include <iostream>
#include <cstdlib>
using namespace std;

int comparisons = 0; 

void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

int partition(int arr[], int low, int high) {
    int pivot = arr[low]; // Use the first element as the pivot
    int start = low;
    int end = high;

    while (start < end) {
        // Increment start while arr[start] <= pivot
        while (arr[start] <= pivot && start < high) {
            comparisons++;
            start++;
        }

        // Decrement end while arr[end] > pivot
        while (arr[end] > pivot && end > low) {
            comparisons++;
            end--;
        }

        // Swap elements if start < end
        if (start < end) {
            swap(&arr[start], &arr[end]);
        }
    }

    // Swap the pivot with arr[end]
    swap(&arr[low], &arr[end]);
    return end; // Return the pivot's final position
}


int partition_r(int arr[], int low, int high) {
    int random = low + rand() % (high - low + 1); 
    swap(&arr[random], &arr[high]);
    return partition(arr, low, high);
}

void quickSort(int arr[], int low, int high) {
    if (low < high) {
        int pi = partition_r(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

void printArray(int arr[], int size) {
    for (int i = 0; i < size; i++)
        cout << arr[i] << " ";
    cout << endl;
}

int main() {
    int n;
    cout << "Enter the size of the array: ";
    cin >> n;
    int arr[n];
    cout << "Enter the elements of the array: ";
    for (int i = 0; i < n; i++) {
        cin >> arr[i];
    }

    quickSort(arr, 0, n - 1);

    cout << "Sorted array: ";
    printArray(arr, n);

    cout << "Number of comparisons: " << comparisons << endl;

    return 0;
}


