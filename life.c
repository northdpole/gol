#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
 * The game of life
 * input:
 * 	argv[1]: Width(N)
 * 	argv[2]: Height(M)
 * 	argv[3]: input filename
 *
 *  file has the following structure:
 *  x-coordinates of alive cell, y-coordinates of alive cell\n
 * e.g.
 * 2,5
 * 3,4
 * 3,3
 * 1,2
 * 150000,23452345
 */

int main(int argc,char **argv){

	char *bfile,**board,*buff,*substr;
	char d;
	int i,j,x,y,n,m;
	FILE *input;

	if( argc <= 3 ){
		printf("Usage: executable *Board Width* *Board Height* *Input Filename*\n");
		exit(-1);
	}

	n = atoi( argv[1] );
	m = atoi( argv[2] );
	bfile = argv[3];

	//printf("%s\n",bfile);

	*board = malloc(m*sizeof( char* ));
	if(*board == NULL){
		printf("Malloc error\n");
		exit(-1);
	}
	for(i = 0; i < m; i++ ){
		board[i] = malloc(n*sizeof(char));
		if(board[i] == NULL){
			printf("Malloc error\n");
			exit(-1);
		}
		memset(board[i], '-', n);
	}
	input = fopen(bfile,"r");
	 if (input == NULL)
		printf("File doesn't exist\n");
	buff = malloc(47*sizeof(char));
	if (buff == NULL){
		printf("Malloc error\n");
		goto exit;
	}

	while (fgets (buff, 47, input) != NULL){
		substr = strtok(buff,",");
	//puts(substr);
		x = atoi(substr);
		if (x > m )
			goto exit;
		substr = strtok(NULL,"\n");
	//puts(substr);
		y = atoi(substr);
		if (y > n )
			goto exit;
		board[x][y] = 'A';
	//printf("%d %d \n", x, y);
	}

	if( m < 40 && n < 80){
		printf("Initial Board: \n");
		for (i = 0; i < m; i++){
			for (j = 0; j < n; j++){
				printf("%c ", board[i][j]);
			}
			printf("\n");
		}
	}
	/*Play life*/
exit:
	for(i = 0; i < m; i++ ){
		free(board[i]);
	}
	free(buff);
	//free(*board);
	printf("Out\n");
return 0;
}
