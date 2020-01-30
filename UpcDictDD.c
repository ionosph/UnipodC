/*
 *  UpcDictDD.c
 *  UnipodC
 *
 *  Created by ionosph on 12/01/27.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcDictDD.h"

#define GETCOLOR(node) ((node)?((node)->color):BLACK)

static const int BLACK = 0;
static const int RED = 1;


static UpcRBNodeDD *UpcRBNodeDD_init(keyType key, valType val, int color)
{
    UpcRBNodeDD *self = (UpcRBNodeDD *)malloc(sizeof(UpcRBNodeDD));

    if (self) {
        self->key = key;
        self->val = val;
        self->color = color;
        self->left = NULL;
        self->right = NULL;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcRBNodeDD_init");
    }
    return self;
}

static void UpcRBNodeDD_dealloc(UpcRBNodeDD *self)
{
    if (self) {
        UpcRBNodeDD_dealloc(self->left);
        UpcRBNodeDD_dealloc(self->right);
        free(self);
    }
}

static UpcRBNodeDD *rotateRight(UpcRBNodeDD *node)
{
    UpcRBNodeDD *lnode = node->left;

    node->left = lnode->right;
    lnode->right = node;
    lnode->color = node->color;
    node->color = RED;
    return lnode;
}

static UpcRBNodeDD *rotateLeft(UpcRBNodeDD *node)
{
    UpcRBNodeDD *rnode = node->right;

    node->right = rnode->left;
    rnode->left = node;
    rnode->color = node->color;
    node->color = RED;
    return rnode;
}

/*
static Bool searchNode(UpcRBNodeDD *node, keyType key, valType *val)
{
    while (node) {
        if (key == node->key) {
            *val = node->val;
            return TRUE;
        }
        else if (key < node->key) {
            node = node->left;
        }
        else {
            node = node->right;
        }
    }
    return FALSE;
}
*/

static void splitNode(UpcRBNodeDD *node)
{
    node->color = RED;
    node->left->color = BLACK;
    node->right->color = BLACK;
}

static UpcRBNodeDD *balanceInsertLeft(UpcRBNodeDD *node, Bool *flg)
{
    if (*flg) {
        return node;
    }
    if (node->color == BLACK) {
        *flg = TRUE;
        if (GETCOLOR(node->left->right) == RED) {
            node->left = rotateLeft(node->left);
        }
        if (GETCOLOR(node->left->left) == RED) {
            if (GETCOLOR(node->right) == RED) {
                splitNode(node);
                *flg = FALSE;
            }
            else {
                node = rotateRight(node);
            }
        }
    }
    return node;
}

static UpcRBNodeDD *balanceInsertRight(UpcRBNodeDD *node, Bool *flg)
{
    if (*flg) {
        return node;
    }
    if (node->color == BLACK) {
        *flg = TRUE;
        if (GETCOLOR(node->right->left) == RED) {
            node->right = rotateRight(node->right);
        }
        if (GETCOLOR(node->right->right) == RED) {
            if (GETCOLOR(node->left) == RED) {
                splitNode(node);
                *flg = FALSE;
            }
            else {
                node = rotateLeft(node);
            }
        }
    }
    return node;
}

static UpcRBNodeDD *insertNode(UpcRBNodeDD *node, keyType key, valType val, Bool *flg)
{

    if (!node) {
        *flg = FALSE;
        return UpcRBNodeDD_init(key, val, RED);
    }
    if (key < node->key) {
        node->left = insertNode(node->left, key, val, flg);
        return balanceInsertLeft(node, flg);
    }
    else if (key > node->key) {
        node->right = insertNode(node->right, key, val, flg);
        return balanceInsertRight(node, flg);
    }
    node->val = val;
    *flg = TRUE;
    return node;
}


static UpcRBNodeDD *balanceRight(UpcRBNodeDD *node, Bool *flg)
{
    if (*flg) {
        return node;
    }
    if (GETCOLOR(node->left->left) == BLACK && GETCOLOR(node->left->right) == BLACK) {
        if (node->left->color == BLACK) {
            node->left->color = RED;
            if (node->color == BLACK) {
                *flg = FALSE;
                return node;
            }
            node->color = BLACK;
        }
        else {
            node = rotateRight(node);
            *flg = FALSE;
            node->right = balanceRight(node->right, flg);
        }
    }
    else {
        if (GETCOLOR(node->left->right) == RED) {
            node->left = rotateLeft(node->left);
        }
        node = rotateRight(node);
        node->right->color = BLACK;
        node->left->color = BLACK;
    }
    *flg = TRUE;
    return node;
}

static UpcRBNodeDD *balanceLeft(UpcRBNodeDD *node, Bool *flg)
{

    if (*flg) {
        return node;
    }
    if (GETCOLOR(node->right->left) == BLACK && GETCOLOR(node->right->right) == BLACK) {
        if (node->right->color == BLACK) {
            node->right->color = RED;
            if (node->color == BLACK) {
                *flg = FALSE;
                return node;
            }
            node->color = BLACK;
        }
        else {
            node = rotateLeft(node);
            *flg = FALSE;
            node->left = balanceLeft(node->left, flg);
        }
    }
    else {
        if (GETCOLOR(node->right->left) == RED) {
            node->right = rotateRight(node->right);
        }
        node = rotateLeft(node);
        node->left->color = BLACK;
        node->right->color = BLACK;
    }
    *flg = TRUE;
    return node;
}

static void searchMin(UpcRBNodeDD *node, keyType *key, valType *val)
{
    while (node->left)
        node = node->left;
    *key = node->key;
    *val = node->val;
}

static UpcRBNodeDD *deleteNode(UpcRBNodeDD *node, keyType key, Bool *flg)
{
    if (!node) {
        *flg = TRUE;
        return node;
    }
    if (key == node->key) {
        if (!(node->left) && !(node->right)) {
            *flg = (node->color == RED) ? TRUE : FALSE;
            UpcRBNodeDD_dealloc(node);
            return NULL;
        }
        else if (!(node->right)) {
            UpcRBNodeDD *lnode = node->left;

            node->left = NULL;
            UpcRBNodeDD_dealloc(node);
            lnode->color = BLACK;
            *flg = TRUE;
            return lnode;
        }
        else if (!(node->left)) {
            UpcRBNodeDD *rnode = node->right;

            node->right = NULL;
            UpcRBNodeDD_dealloc(node);
            rnode->color = BLACK;
            *flg = TRUE;
            return rnode;
        }
        else {
            keyType k;
            valType v;

            searchMin(node->right, &k, &v);
            node->key = k;
            node->val = v;
            node->right = deleteNode(node->right, k, flg);
            return balanceRight(node, flg);
        }
    }
    else if (key < node->key) {
        node->left = deleteNode(node->left, key, flg);
        return balanceLeft(node, flg);
    }
    else {
        node->right = deleteNode(node->right, key, flg);
        return balanceRight(node, flg);
    }
}

static int checkRBTree(UpcRBNodeDD *node)
{
    int a, b;

    if (node) {
        if (node->color == RED) {
            if (GETCOLOR(node->left) == RED || GETCOLOR(node->right) == RED) {
                fprintf(stderr, "RBTree Check Error 1.\n");
                exit(1);
            }
        }
        a = checkRBTree(node->left);
        b = checkRBTree(node->right);
        if (a != b) {
            fprintf(stderr, "RBTree Check Error 2.\n");
            exit(1);
        }
        if (node->color == BLACK) {
            a++;
        }
        return a;
    }
    return 0;
}

static void printRBTree(UpcRBNodeDD *node, int n)
{
    int i;

    if (node) {
        printRBTree(node->left, n + 1);
        for (i = 0; i < n; i++)
            printf("    ");
        printf("%c(%10.3e:%10.3e)\n", ((node->color)?'R':'B'), node->key, node->val);
        printRBTree(node->right, n + 1);
    }
}

Bool test_UpcRBNodeDD(void)
{
    UpcRBNodeDD *root = NULL;
    int i, n = 8;
    Bool flg;

    for (i = 0; i < n; i++) {
        printf("-------- insert %d\n", i);
        root = insertNode(root, i, 0.0, &flg);
        root->color = BLACK;
        printRBTree(root, 0);
        checkRBTree(root);
    }
    return TRUE;
}

//----------------------------------------

static void UpcDictDD_dealloc(UpcDictDD *self)
{
    UpcRBNodeDD_dealloc(self->root);
    free(self);
}

UpcDictDD *UpcDictDD_init(void)
{
    UpcDictDD *self = (UpcDictDD *)malloc(sizeof(UpcDictDD));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcDictDD_dealloc;
        self->root = NULL;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcDictDD_init");
    }
    return self;
}


Bool UpcDictDD_search(UpcDictDD *self, keyType key, valType *val)
{
    UpcRBNodeDD *node = self->root;

    while (node) {
        if (key == node->key) {
            *val = node->val;
            return TRUE;
        }
        else if (key < node->key) {
            node = node->left;
        }
        else {
            node = node->right;
        }
    }
    return FALSE;
}

void UpcDictDD_add(UpcDictDD *self, keyType key, valType val)
{
    Bool flg;

    self->root = insertNode(self->root, key, val, &flg);
    self->root->color = BLACK;
}

void UpcDictDD_delete(UpcDictDD *self, keyType key)
{
    Bool flg;

    deleteNode(self->root, key, &flg);
}


static void printUpcDictDDNode(UpcRBNodeDD *node, FILE *ostream)
{
    if (node) {
        printUpcDictDDNode(node->left, ostream);
        fprintf(ostream, "  %20.12e:%20.12e\n", node->key, node->val);        
        printUpcDictDDNode(node->right, ostream);
    }
}

void UpcDictDD_print(UpcDictDD *self, FILE *ostream)
{
    fprintf(ostream, "     key             ,    val\n");
    printUpcDictDDNode(self->root, ostream);
}

Bool TestUpcDictDD(void)
{
    UpcDictDD *t = UpcDictDD_init();
    int i;
    double x;

    for (i = 0; i < 10; i++) {
        UpcDictDD_add(t, i, i * i);
    }
    for (i = 1; i >= -10; i--) {
        UpcDictDD_add(t, i, i * i);
    }
    for (x = 0.5; x < 20; x += 2) {
        UpcDictDD_add(t, x, x * x);
    }
    checkRBTree(t->root);
    UpcDictDD_print(t, stdout);
    UpcDictDD_add(t, 5.0, 555.0);
    checkRBTree(t->root);
    UpcDictDD_print(t, stdout);
    UpcDictDD_add(t, 0.01, 222.0);
    checkRBTree(t->root);
    UpcDictDD_print(t, stdout);
    printf("%s\n", (UpcDictDD_search(t, 6.0, &x) ? "TRUE" : "FALSE"));
    printf("x = %.12e\n", x);
    printRBTree(t->root, 0);

    for (i = 0; i < 10; i++) {
        UpcDictDD_delete(t, i);
    }
    checkRBTree(t->root);
    printf("-----\n");
    UpcDictDD_print(t, stdout);
    //printRBTree(t->root, 0);

    UpcDictDD_delete(t, 0.123);
    //printRBTree(t->root, 0);

    UpcObj_RELEASE(t);
    return TRUE;
}




