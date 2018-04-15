//  NO LONGER standalone code
/*
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
*/


/*
 * This file is to get together all the helper function for the exact meniscus curve from the perspective of IO (NON-NUMERICAL)
 * To be clear, there are some helper function in imkcurve.c. That's more of numerical side. (NUMERICAL)
 * API-like, but not entirely
 */

#include"iprotos.h"
struct DataItem* hashArray[MSIZE];
struct DataItem* dummyItem;
struct DataItem* item;

EXPORT int hashCode(int key)
{
    return key;
}

EXPORT struct DataItem *search(int key)
{
    //get the hash
    int hashIndex = hashCode(key);

    //move in array until an empty
    while(hashArray[hashIndex] != NULL) {

        if(hashArray[hashIndex]->key == key)
            return hashArray[hashIndex];

        //go to next cell
        ++hashIndex;
    }
    return NULL;
}

EXPORT void insert(int key,double data)
{
    struct DataItem *item = (struct DataItem*) malloc(sizeof(struct DataItem));
    item->data = data;
    item->key = key;

    //get the hash
    int hashIndex = hashCode(key);

    //move in array until an empty or deleted cell
    while(hashArray[hashIndex] != NULL && hashArray[hashIndex]->key != -1) {
        //go to next cell
        ++hashIndex;

    }

    hashArray[hashIndex] = item;
    //printf("in func %s, hashIndex = %d data = %lf key = %d\n", __func__, hashIndex, data, key); // TODO ** : Print statement could be commented out or leave it in debugging mode. This is a minor deal
}

struct DataItem* delete(struct DataItem* item)
{
    int key = item->key;

    //get the hash
    int hashIndex = hashCode(key);

    //move in array until an empty
    while(hashArray[hashIndex] != NULL) {

        if(hashArray[hashIndex]->key == key) {
            struct DataItem* temp = hashArray[hashIndex];

            //assign a dummy item at deleted position
            hashArray[hashIndex] = dummyItem;
            return temp;
        }

        //go to next cell
        ++hashIndex;
    }

    return NULL;
}

// this function serves as a double check when data was read and passed into Hash Table
// label: debugging
EXPORT void display()
{
    int i = 0;

    // TODO: iterator upper limit was changed to MN instead of MSIZE due to the fact that we are currently concerning MN points ONLY.
    for(i = 0; i<MN; i++) {

        if(hashArray[i] != NULL)
            printf(" (%d,%lf)",hashArray[i]->key,hashArray[i]->data);
        else
            printf(" ~~ ");
    }

    printf("\n");
}

// this function should be called once ONLY for the purpose of reading exp folder
//
EXPORT void readMeniscusFromTxtFile(int CA, int EXP) // take contact angle and its experiment number as input
{
    printf("%s started reading process\n", __func__);
    dummyItem = (struct DataItem*) malloc(sizeof(struct DataItem));
    dummyItem->data = -1.;
    dummyItem->key = -1;


    // READ input from FILE
    FILE *fp;
    char filename[20];
    char dirname[20];
    //int CA = 20;  // contact angle arg
    //int EXP = 105; // experiment arg
    double x, z;
    double L = 0; // Lower boundary
    double U = 1.875; // Upper boundary
    double spacing = (U-L)/ (MN-1); // This is the spacing used in iFluid simulation code

    printf("spacing is %lf\n", spacing);
    sprintf(dirname, "sim%d", EXP);
    sprintf(filename, "exp%d-%d.txt", EXP, CA);
    char *newf = strcat(strcat(dirname, "/"), filename);
    fp = fopen(newf,"r");
    for (int i = 0; i < MN; i++)
    {
        fscanf(fp, "%lf %lf\n", &x, &z);
        //printf("%lf %lf \n", x, z); // passed test and commented out
        int k = (int)(round((x - L) / spacing)); // FIXME: round is a good choice
        //printf("k = %d x = %lf L = %lf, ceil_ratio = %lf ratio = %lf\n", k, x, L, ceil((x-L) / spacing), (x-L)/spacing); // passed test and commented out
        insert(k, z);
    }
    fclose(fp);
    display();
    printf("Meniscus was read during initial setup\n");
}

// this function is to match up the coordinate produced in FronTier with exp folders which is from OCTAVE solver
// L: is the leftmost domain limit
// spacing: is the mesh grid size
// q: is the x coordinate
// function type: float, return the coordinate of adjusted z_coordinate
EXPORT boolean checkNumericalMeniscusforZ(double q, double L, double spacing)
{
    //display(); // This is a test for the purpose of that hash table is in memory. // TODO: should be commented out when testing is finished
    // convert coordinate to key of hash table
    int q_index = (int)(round((q-L)/ spacing)); // TODO: still testing here

    // hash table lookup here:
    if (search(q_index))
    {
        //printf("%lf was matched up!\n", q);
        return YES;
    }
    else
    {
        //printf("No coordinates were matched for %lf!\n", q);
        return NO;
    }
}
