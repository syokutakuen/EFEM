//
//  雑役関数群
//

//  エラー報告関数
void error (int code, char *msg);

//	番号から配列上の位置を探す関数群
int search_einfo (int n);
int search_part (int n);
int search_node (int n);
int search_elem (int n);
int search_mat (int n);

//	リストに当該データーを追加
void attach_node (int *ip, int n);
void attach_elem (int *ip, int n);
void attach_part (int *ip, int n);

