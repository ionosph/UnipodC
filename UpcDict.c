/*
 *  UpcDict.c
 *  UnipodC
 *
 *  Created by ionosph on 12/01/27.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "UpcDict.h"

#define GETCOLOR(node) ((node)?((node)->color):BLACK)

static const int BLACK = 0;
static const int RED = 1;

static UpcRBNode *UpcRBNode_init(void *key, void *val, int color)
{
    UpcRBNode *self = (UpcRBNode *)malloc(sizeof(UpcRBNode));

    if (self) {
        self->key = key; // (注) 所有権が移る.
        self->val = val; // (注) 所有権が移る.
        self->color = color;
        self->left  = NULL;
        self->right = NULL;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcRBNode_init");
    }
    return self;
}

static void UpcRBNode_dealloc(UpcRBNode *self, void (*delKey)(void *), void (*delVal)(void *))
{
    if (self) {
        UpcRBNode_dealloc(self->left , delKey, delVal);
        UpcRBNode_dealloc(self->right, delKey, delVal);
        delKey(self->key);
        delVal(self->val);
        free(self);
    }
}

static UpcRBNode *rotateRight(UpcRBNode *node)
{
    UpcRBNode *lnode = node->left;

    node->left = lnode->right;
    lnode->right = node;
    lnode->color = node->color;
    node->color = RED;
    return lnode;
}

static UpcRBNode *rotateLeft(UpcRBNode *node)
{
    UpcRBNode *rnode = node->right;

    node->right = rnode->left;
    rnode->left = node;
    rnode->color = node->color;
    node->color = RED;
    return rnode;
}

/*
static Bool searchNode(UpcRBNode *node, int (*compKey)(const void *, const void *), void *key, void **val)
{
    int cmp;

    while (node) {
        cmp = compKey(key, node->key);
        if (cmp == 0) {
            *val = node->val;
            return TRUE;
        }
        else if (cmp < 0) {
            node = node->left;
        }
        else {
            node = node->right;
        }
    }
    return FALSE;
}
*/

static void splitNode(UpcRBNode *node)
{
    node->color = RED;
    node->left->color = BLACK;
    node->right->color = BLACK;
}

static UpcRBNode *balanceInsertLeft(UpcRBNode *node, Bool *flg)
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

static UpcRBNode *balanceInsertRight(UpcRBNode *node, Bool *flg)
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

static UpcRBNode *insertNode(UpcRBNode *node, int (*compKey)(const void *, const void *), void (*delKey)(void *), void (*delVal)(void *), void *key, void *val, Bool *flg)
{
    int cmp;

    if (!node) {
        *flg = FALSE;
        return UpcRBNode_init(key, val, RED);
    }
    cmp = compKey(key, node->key);
    if (cmp < 0) {
        node->left = insertNode(node->left, compKey, delKey, delVal, key, val, flg);
        return balanceInsertLeft(node, flg);
    }
    else if (cmp > 0) {
        node->right = insertNode(node->right, compKey, delKey, delVal, key, val, flg);
        return balanceInsertRight(node, flg);
    }
    // 重複したキーがあったとき
    delKey(key);
    delVal(node->val);
    node->val = val;
    *flg = TRUE;
    return node;
}

static UpcRBNode *balanceRight(UpcRBNode *node, Bool *flg)
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

static UpcRBNode *balanceLeft(UpcRBNode *node, Bool *flg)
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

static UpcRBNode *searchMin(UpcRBNode *node)
{
    while (node->left)
        node = node->left;
    return node;
}

static UpcRBNode *deleteNode(UpcRBNode *node, int (*compKey)(const void *, const void *), void (*delKey)(void *), void (*delVal)(void *), void *key, Bool *flg)
{
    int cmp;

    if (!node) {
        *flg = TRUE;
        return node;
    }
    cmp = compKey(key, node->key);
    if (cmp == 0) {
        if (!(node->left) && !(node->right)) {
            *flg = (node->color == RED) ? TRUE : FALSE;
            UpcRBNode_dealloc(node, delKey, delVal);
            return NULL;
        }
        else if (!(node->right)) {
            UpcRBNode *lnode = node->left;

            node->left = NULL;
            UpcRBNode_dealloc(node, delKey, delVal);
            lnode->color = BLACK;
            *flg = TRUE;
            return lnode;
        }
        else if (!(node->left)) {
            UpcRBNode *rnode = node->right;

            node->right = NULL;
            UpcRBNode_dealloc(node, delKey, delVal);
            rnode->color = BLACK;
            *flg = TRUE;
            return rnode;
        }
        else {
            UpcRBNode *t = searchMin(node->right);
            void *k = t->key;
            void *v = t->val;

            t->key = node->key;
            t->val = node->val;
            node->key = k;
            node->val = v;
            node->right = deleteNode(node->right, compKey, delKey, delVal, key, flg);
            return balanceRight(node, flg);
        }
    }
    else if (cmp < 0) {
        node->left = deleteNode(node->left, compKey, delKey, delVal, key, flg);
        return balanceLeft(node, flg);
    }
    else {
        node->right = deleteNode(node->right, compKey, delKey, delVal, key, flg);
        return balanceRight(node, flg);
    }
}

static int UpcRBNode_count(UpcRBNode *node)
{
    int l, r;

    if (!node)
        return 0;

    l = UpcRBNode_count(node->left);
    r = UpcRBNode_count(node->right);
    return l + r + 1;
}

static int UpcRBNode_serializeKey(UpcRBNode *node, int current, void **keys)
{
    int i;

    if (!node)
        return current;

    i = UpcRBNode_serializeKey(node->left, current, keys);
    keys[i] = node->key;
    return UpcRBNode_serializeKey(node->right, i + 1, keys);
}

static int UpcRBNode_serializeVal(UpcRBNode *node, int current, void **vals)
{
    int i;
    
    if (!node)
        return current;
    
    i = UpcRBNode_serializeVal(node->left, current, vals);
    vals[i] = node->val;
    return UpcRBNode_serializeVal(node->right, i + 1, vals);
}

/*
static int checkRBTree(UpcRBNode *node)
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

static void printRBTree(UpcRBNode *node, int n)
{
    int i;

    if (node) {
        printRBTree(node->left, n + 1);
        for (i = 0; i < n; i++)
            printf("    ");
        printf("%c(%p:%p)\n", ((node->color)?'R':'B'), node->key, node->val);
        printRBTree(node->right, n + 1);
    }
}
*/

//----------------------------------------

static void UpcDict_dealloc(UpcDict *self)
{
    UpcRBNode_dealloc(self->root, self->delKey, self->delVal);
    free(self);
}

UpcDict *UpcDict_init(int (*compKey)(const void *, const void *), void (*delKey)(void *), void (*delVal)(void *))
{
    UpcDict *self = (UpcDict *)malloc(sizeof(UpcDict));

    if (self) {
        self->_refCount = 1;
        self->_dealloc = UpcDict_dealloc;
        self->root = NULL;
        self->compKey = compKey;
        self->delKey = delKey;
        self->delVal = delVal;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcDict_init");
    }
    return self;
}


Bool UpcDict_search(UpcDict *self, const void *key, void **val)
{
    UpcRBNode *node = self->root;
    int cmp;

    while (node) {
        cmp = self->compKey(key, node->key);
        if (cmp == 0) {
            *val = node->val;
            return TRUE;
        }
        else if (cmp < 0) {
            node = node->left;
        }
        else {
            node = node->right;
        }
    }
    return FALSE;
}

/* この関数に渡したkey, valは所有権がUpcDictオブジェクトに移る.
 * 必要なくなったときの解放などはUpcDictオブジェクトが受け持つ.
 * UpcDictオブジェクトに解放されては困る場合はRETAINしておくこと.
 */
void UpcDict_add(UpcDict *self, void *key, void *val)
{
    Bool flg;

    self->root = insertNode(self->root, self->compKey, self->delKey, self->delVal, key, val, &flg);
    self->root->color = BLACK;
}

void UpcDict_delete(UpcDict *self, void *key)
{
    Bool flg;

    self->root = deleteNode(self->root, self->compKey, self->delKey, self->delVal, key, &flg);
    self->root->color = BLACK;
}

static void printUpcDictNode(UpcRBNode *node, FILE *ostream)
{
    if (node) {
        printUpcDictNode(node->left, ostream);
        fprintf(ostream, "  %s: %p\n", node->key, node->val);        
        printUpcDictNode(node->right, ostream);
    }
}

void UpcDict_print(UpcDict *self, FILE *ostream)
{
    fprintf(ostream, "     key             ,    val\n");
    printUpcDictNode(self->root, ostream);
}

/* keyの配列を返す.
 * この関数が返した配列は, 使い終わったら解放すること.
 */
void **UpcDict_keys(UpcDict *self)
{
    int n = UpcRBNode_count(self->root);
    void **keys = (void **)malloc(sizeof(void *) * (n + 1));

    if (keys) {
        UpcRBNode_serializeKey(self->root, 0, keys);
        keys[n] = NULL;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcDict_keys");
    }
    return keys;
}

/* valの配列を返す.
 * この関数が返した配列は, 使い終わったら解放すること.
 */
void **UpcDict_vals(UpcDict *self)
{
    int n = UpcRBNode_count(self->root);
    void **vals = (void **)malloc(sizeof(void *) * (n + 1));

    if (vals) {
        UpcRBNode_serializeVal(self->root, 0, vals);
        vals[n] = NULL;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcDict_vals");
    }
    return vals;
}

// この関数が返した配列は, 使い終わったら解放すること.
/*
void **UpcDict_items(UpcDict *self)
{
    int n = UpcRBNode_count(self->root);
    void **items = (void **)malloc(sizeof(void *) * (2 * n + 1));

    if (keys) {
        UpcRBNode_serializeItems(self->root, 0, items);
        items[2 * n] = NULL;
    }
    else {
        UpcERRHANDLER(UpcE_MemoryError, "in UpcDict_keys");
    }
    return items;
}
*/

