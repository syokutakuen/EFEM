//	雑役関数群

#include <stdlib.h>
#include "typedef.h"
#include "memfree.h"

//	エラー報告関数
//	エラー報告後、確保したメモリを開放する。
void error (int code, char *msg)
{
	fprintf (stderr, "%s at %d\n", msg, nline);
	free_data ();
	exit (code);
}

//	番号から配列上の位置を探す関数群
int search_einfo (int n)
{
	int i;
	for (i = 0; i < neinfo; i++)
		if (einfo [i].number == n) break;
	return i;
}

int search_part (int n)
{
	int i;
	for (i = 0; i < npart; i++)
		if (part [i].number == n) break;
	return i;
}

int search_node (int n)
{
	int i;
	for (i = 0; i < ntnode; i++)
		if (node [i].number == n) break;
	return i;
}

int search_elem (int n)
{
	int i;
	for (i = 0; i < ntelem; i++)
		if (elem [i].number == n) break;
	return i;
}

int search_mat (int n)
{
	int i;
	for (i = 0; i < nmat; i++)
		if (mat [i]->number == i) break;
	return i;
}

//	リストに当該データーを追加
void attach_node (int *ip, int n)
{
	int ptr = search_node (n);
	if (ptr == ntnode) error (-99, "Illegal Node number");
	*ip = ptr;
}

void attach_elem (int *ip, int n)
{
	int ptr = search_elem (n);
	if (ptr == ntelem) error (-99, "Illegal elem number");
	*ip = ptr;
}

void attach_part (int *ip, int n)
{
	int ptr = search_part (n);
	if (ptr == npart) error (-99, "Illegal part number");
	*ip = ptr;
}


