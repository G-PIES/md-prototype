#define IN
#define OUT

kernel void MatrixAdd( IN global const float *dA, IN global const float *dB, IN global int *dMW, OUT global float *dC )
{
	// [dA] is dMW x dMW
	// [dB] is dMW x dMW
	// [dC] is dMW x dMW
	// but all the matrixs' rows are really linear in memory

	int mw = *dMW;
	int crow = get_global_id( 0 );
	int ccol = get_global_id( 1 );

	/*
	int aindex = crow * mw;			// a[i][0]
	int bindex = ccol;				// b[0][j]
	*/
	int cindex = crow * mw + ccol;	// c[i][j]

	/*
	float cij = 0.;
	for( int k = 0; k < mw; k++ )
	{
		cij += dA[aindex] + dB[bindex];
		aindex++;
		bindex += mw;
	}
	dC[cindex] = cij;
	*/
	dC[cindex] = dA[cindex] + dB[cindex];
}