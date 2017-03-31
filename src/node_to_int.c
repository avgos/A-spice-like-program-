#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "node_to_int.h"
#define TABLE_SIZE 104179 //prime
struct node * table[TABLE_SIZE];
char ** _node_names;
int _node_names_count;
int _node_names_size;

#define LOWERCASE 1


uint32_t jenkins_one_at_a_time_hash(char *key, size_t len)
{
    uint32_t hash, i;
    for(hash = i = 0; i < len; ++i)
    {
        #if LOWERCASE
            hash += tolower(key[i]);
        #else
            hash += key[i];
        #endif
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

char * get_nodename(int index){
    if(index < _node_names_size)
        return _node_names[index];
    else
        return NULL;
}
void _add_name(char * name){
    if(_node_names_count + 1 == _node_names_size){
        int tmp = _node_names_size;
        _node_names_size += 256*8;
        _node_names = realloc(_node_names, sizeof(char*)*_node_names_size);
        for(int i=tmp; i<_node_names_size;i++)
            _node_names[i] = NULL;
    }
    _node_names[_node_names_count] = strdup(name);
    _node_names_count++;
}


void init_hash(){
    _node_names_count = 0;
    _node_names_size = 512;
    //_node_names = malloc(sizeof(char *)*_node_names_size);
    _node_names = calloc(sizeof(char *), _node_names_size);

    table_counter = 1;
    int i = 0;
    for(i=0;i<TABLE_SIZE;i++)
        table[i] = NULL;
}


int get_key(char * val){
    struct node * temp;
    uint32_t key;
    temp = lookup(val,&key);
    if(temp != NULL){
        return temp->key;
    }
    return -1;
}

struct node * lookup(char * val, uint32_t * key){
    *key = jenkins_one_at_a_time_hash(val,strlen(val));
    uint32_t _key = *key % TABLE_SIZE;
    struct node *temp;
    temp = table[_key];
    while(temp != NULL){
        #if LOWERCASE
            if(strcasecmp(temp->val,val) == 0)
        #else
            if(strcmp(temp->val,val) == 0)
        #endif
            return temp;
        temp = temp->next;
    }
    return NULL;
}

uint32_t add_node(char * val){
    static struct node * next_key = NULL;
    //return 1;
    if(val[0] == '0' && val[1] == '\0')
        return 0;
    struct node * temp;
    uint32_t key;
    temp = lookup(val,&key);
    if(temp != NULL){
        return temp->key;
    }
    //printf("ADD NODE: %s %d\n",val,table_counter );
    _add_name(val);
    key = key % TABLE_SIZE;
    temp = malloc(sizeof(struct node));
    if(!temp){//error
        fprintf(stderr, "MALLOC FAILED FILE: %s FUNCTION: %s LINE: %d\n",__FILE__,__FUNCTION__,__LINE__);
        exit(-1);
    }
    temp->val = strdup(val);
    temp->key = table_counter++;
    temp->next = table[key];
    temp->next_key = next_key;
    next_key = temp;
    table[key] = temp;
    return temp->key;
}
