#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
 * The game of life
 * input:
 * 	argv[1]: Width(N)
 * 	argv[2]: Height(M)
 * 	argv[3]: input filename
 */
int main(int argc,char **argv){

	char *bfile,**board;
	char d;
	int i,j,x,y,n,m;

	if( argc <= 3 ){
		printf("Usage: executable *Board Width* *Board Height* *Input Filename*\n");
		exit(-1);
	}

	n = atoi( argv[1] );
	m = atoi( argv[2] );
	bfile = argv[3];

	//printf("%s\n",bfile);

	*board = malloc( m*sizeof( char* ));
	if(*board == NULL)
	exit(-1);
	for(i = 0; i < m; i++ ){
		board[i] = malloc(n*sizeof(char));
		if(board[i] == NULL)
		exit(-1);
	}


}
