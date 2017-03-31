#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

struct node{
    char *val;
    double key;
    struct node * next;
};


//#define THRESHOLD 0.0000000000000001
#define THRESHOLD 0.01

#define TABLE_SIZE 104179 //prime
struct node * table[TABLE_SIZE];

uint32_t add_node(char * val, char *valt, double key);
struct node * lookup(char * val, int len, uint32_t * key);
void init_table();

unsigned int open_files(int *_fd1, int *_fd2, char **_map1, int * map1_size, char **_map2, int * map2_size, const char *filename_1, const char * filename_2);

int main(int argc, const char *argv[])
{    
    if(argc != 3){
        fprintf(stderr, "USAGE: %s input1 input2\n",argv[0] );
        return -1;
    }

    init_table();
    int fd1, fd2;
    int size_1, size_2;
    char *map1, *map2;
    open_files(&fd1, &fd2, &map1, &size_1, &map2, &size_2, argv[1], argv[2]);

    int start_word = 0;
    int end_word = 0;
    int i = 0;
    int counter = 0;
    unsigned int total_mem = 0;
    int size;
    size = size_1;
    char * map;
    map = map1;
    while(i< size){
        while(map[i] == ' '  && i < size) i++;
        while(map[i] == '\n' && i < size) i++;
        start_word = i;
        end_word = i;
        while(map[end_word] != ' ' && end_word < size) end_word++;
    
        if(end_word == i){
            //zero length word.. go to the end of line
            while(map[i] == '\n' && i < size) i++;
         }else{
            i = end_word;

        }
        while(map[i] == ' ') i++;
        int start_num = i;
        int end_num = i;
        while(map[end_num] != '\n' && end_num < size) end_num++;
        if(end_num == i){
            //error 
            while(map[i] == '\n' && i < size) i++;
        }else{
            char * end = &map[end_num];
            total_mem += sizeof(struct node) + end_word - start_word;
            add_node(&map[start_word], &map[end_word], strtod(&map[start_num], &end));
            counter ++;
            i = end_num;
        }
    }

    i = 0;
    size = size_2;
    map = map2;
    int threshold_violations = 0;
    double max_norm =  -2.225E-307;
    double relative_error = 0.0;
    double norm_relative_error = 0.0;

    int failed_nodes = 0;
    int good_nodes = 0;
    int total_nodes = 0;
    while(i< size){

        while(map[i] == ' '  && i < size) i++;
        while(map[i] == '\n' && i < size) i++;
        start_word = i;
        end_word = i;
        while(map[end_word] != ' ' && end_word < size) end_word++;
    
        if(end_word == i){
            //zero length word.. go to the end of line
            while(map[i] == '\n' && i < size) i++;
         }else{
            i = end_word;

        }
        while(map[i] == ' ') i++;
        int start_num = i;
        int end_num = i;
        while(map[end_num] != '\n' && end_num < size) end_num++;
        if(end_num == i){
            //error 
            while(map[i] == '\n' && i < size) i++;
        }else{
            char * end = &map[end_num];
            total_mem += sizeof(struct node) + end_word - start_word;
            uint32_t key;
            total_nodes++;
            struct node * ht_node =  lookup(&map[start_word], end_word - start_word, &key);
            if(ht_node == NULL){
                failed_nodes++;
            }else{
                good_nodes++;
                if(ht_node->key != 0.0){
                    relative_error += fabs( 1 - (strtod(&map[start_num], &end)/ht_node->key) );
                    norm_relative_error += (1 - (strtod(&map[start_num], &end)/ht_node->key))*(1 - (strtod(&map[start_num], &end)/ht_node->key));
                }
                if(fabs(ht_node->key - strtod(&map[start_num], &end)) > max_norm)
                    max_norm = fabs(ht_node->key - strtod(&map[start_num], &end));
                threshold_violations += fabs(ht_node->key - strtod(&map[start_num], &end)) > THRESHOLD;
            }
            counter ++;
            i = end_num;
        }
    }
    norm_relative_error = sqrt(norm_relative_error);
    printf("failed_nodes: %d good_nodes: %d total_nodes: %d\n", failed_nodes, good_nodes, total_nodes);
    printf("max_norm: %.20e\n",max_norm);
    printf("relative_error: %.20e %f\n",relative_error,relative_error);
    printf("norm_relative_error: %.20e %f\n",norm_relative_error,norm_relative_error);
    printf("threshold_violations: %d\n", threshold_violations );
    
    // Don't forget to free the mmapped memory
    if (munmap(map1, size_1) == -1)
    {
        close(fd1);
        perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    }
    // Don't forget to free the mmapped memory
    if (munmap(map2, size_2) == -1)
    {
        close(fd2);
        perror("Error un-mmapping the file");
        exit(EXIT_FAILURE);
    }
    // Un-mmaping doesn't close the file, so we still need to do that.
    close(fd1);
    close(fd2);
   
    return 0;
}


unsigned int open_files(int *_fd1, int *_fd2, char **_map1, int * map1_size, char **_map2, int * map2_size, const char *filename_1, const char * filename_2){
    
    int fd1 = open(filename_1, O_RDONLY, (mode_t)0600);
    if (fd1 == -1)
    {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }        
    posix_fadvise(fd1, 0, 0, 1);  // FDADVICE_SEQUENTIAL

    
    struct stat fileInfo1 = {0};
    
    if (fstat(fd1, &fileInfo1) == -1)
    {
        perror("Error getting the file size");
        exit(EXIT_FAILURE);
    }
    
    if (fileInfo1.st_size == 0)
    {
        fprintf(stderr, "Error: File is empty, nothing to do\n");
        exit(EXIT_FAILURE);
    }
    *map1_size = fileInfo1.st_size;
    
    
    char *map1 = mmap(0, fileInfo1.st_size, PROT_READ, MAP_SHARED, fd1, 0);
    if (map1 == MAP_FAILED)
    {
        close(fd1);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
    if(madvise(map1, fileInfo1.st_size, MADV_SEQUENTIAL) == -1)
    {
        perror("madvise failed\n");
        exit(EXIT_FAILURE);
    }

    int fd2 = open(filename_2, O_RDONLY, (mode_t)0600);
    if (fd2 == -1)
    {
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }        
    posix_fadvise(fd2, 0, 0, 1);  // FDADVICE_SEQUENTIAL

    struct stat fileInfo2 = {0};
    
    if (fstat(fd2, &fileInfo2) == -1)
    {
        perror("Error getting the file size");
        exit(EXIT_FAILURE);
    }
    
    if (fileInfo2.st_size == 0)
    {
        fprintf(stderr, "Error: File is empty, nothing to do\n");
        exit(EXIT_FAILURE);
    }
    *map2_size = fileInfo2.st_size;
    
//    printf("File size is %ji\n", (intmax_t)fileInfo1.st_size);
    
    char *map2 = mmap(0, fileInfo2.st_size, PROT_READ, MAP_SHARED, fd2, 0);
    if (map2 == MAP_FAILED)
    {
        close(fd2);
        perror("Error mmapping the file");
        exit(EXIT_FAILURE);
    }
    
    if(madvise(map2, fileInfo2.st_size, MADV_SEQUENTIAL) == -1)
    {
        perror("madvise failed\n");
        exit(EXIT_FAILURE);
    }

    *_fd1 = fd1;
    *_fd2 = fd2;
    *_map1 = map1; 
    *_map2 = map2;
    unsigned int min_size = fileInfo2.st_size;
    if(min_size > fileInfo1.st_size)
        min_size = fileInfo1.st_size;
    return min_size;
}


void init_table(){
    int i = 0;
    for(i=0;i<TABLE_SIZE;i++)
        table[i] = NULL;
}

uint32_t jenkins_one_at_a_time_hash(char *key, size_t len)
{
    uint32_t hash, i;
    for(hash = i = 0; i < len; ++i)
    {
        //#if LOWERCASE
          //  hash += tolower(key[i]);
        //#else
            hash += key[i];
        //#endif
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

struct node * lookup(char * val, int len, uint32_t * key){
    *key = jenkins_one_at_a_time_hash(val,len);
    uint32_t _key = *key % TABLE_SIZE;
    struct node *temp;
    temp = table[_key];
    while(temp != NULL){
     //   #if LOWERCASE
        if(strncasecmp(temp->val,val, len) == 0)
       // #else
          //  if(strcmp(temp->val,val) == 0)
        //#endif
            return temp;
        temp = temp->next;
    }
    return NULL;
}

uint32_t add_node(char * val, char *valt, double _key){
    //return 1;
    if(val[0] == '0' && val[1] == '\0')
        return 0;
    struct node * temp;
    uint32_t key;
    temp = lookup(val,valt - val, &key);
    if(temp != NULL){
        return temp->key;
    }
    //printf("ADD NODE: %s %d\n",val,table_counter );

    key = key % TABLE_SIZE;
    temp = malloc(sizeof(struct node));
    if(!temp){//error
        fprintf(stderr, "MALLOC FAILED FILE: %s FUNCTION: %s LINE: %d\n",__FILE__,__FUNCTION__,__LINE__);
        exit(-1);
    }
    char * tmp;
    tmp = malloc((valt - val)*sizeof(char));
    memcpy(tmp, val, valt-val);

    temp->val = tmp;
    temp->key = _key;
    temp->next = table[key];
    table[key] = temp;
    return temp->key;
}
