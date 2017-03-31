#ifndef NODE_TO_INT_H
#define NODE_TO_INT_H
#include <stdint.h>
struct node{
    char *val;
    uint32_t key;
    struct node * next_key;
    struct node * next;
};
uint32_t table_counter; //n
void init_hash();
uint32_t add_node(char * val);
struct node * lookup(char * val, uint32_t * key);
char * get_nodename(int index);
int get_key(char * val);
#endif