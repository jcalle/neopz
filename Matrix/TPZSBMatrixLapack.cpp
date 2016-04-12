/**
 * @file
 * @brief Contains the implementation of the TPZSBMatrixLapack methods.
 */

#include <math.h>
#include <stdlib.h>

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZSBMatrixLapack.h"

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzsbmatrix"));
#endif

using namespace std;

/*******************/
/*** TPZSBMatrixLapack ***/

/**************************** PUBLIC ****************************/

/*****************************/
/*** Construtor (int) ***/

//ok
template<class TVar> TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack() : TPZMatrix<TVar>(0,0) {
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
  fDiag = NULL;
  fBand = 0;
}
//ok
template<class TVar>
TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack( long dim, long band )
: TPZMatrix<TVar>( dim, dim )
{
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
	fBand = ( band > (dim - 1) ? (dim - 1) : band );
	fDiag = new TVar[Size()];
	if ( fDiag == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSBMatrixLapack( dim ) <Error creating Matrix>" );
	
	Zero();
}

//ok
template<class TVar>
TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack(const TPZSBMatrixLapack<TVar> &A ) : TPZMatrix<TVar>(A)  {
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
  Copy(A);
}


/**************/
/*** PutVal ***/
//ok
template <class TVar>
int
TPZSBMatrixLapack<TVar>::PutVal(const long r,const long c,const TVar& value )
{
  // initialises row and col in order to work with upper triangular matrix
  long row(r),col(c);
  if ( row > col )
    this->Swap( &row, &col );
  
  long index;
  if ( (index = col-row) > fBand )
  {
#ifdef PZDEBUG
    if (value != 0. ) {
      DebugStop();
    }
#endif
    return( 0 );        // Element out of bounds
  }
  const long dim = this->fCol;//square matrix
  fDiag[ dim * (fBand - index) + col ] = value;
  this->fDecomposed = 0;
  return( 1 );
}

/**************/
/*** GetVal ***/
//ok
template<class TVar>
const TVar
&TPZSBMatrixLapack<TVar>::GetVal(const long r,const long c ) const
{
	
	// initialises row and col in order to work with upper triangular matrix
  long row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	long index;
	if ( (index = col-row) > fBand )
		return( this->gZero );        // Element out of bounds
	const long dim = this->fCol;//square matrix
	return( fDiag[ dim * (fBand - index) + col ] );
}

/*************/
/*** Print ***/
//ok
template<class TVar>
void
TPZSBMatrixLapack<TVar> ::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const
{
	out.width( 8 );
	out.precision( 4 );
	
	out << "Writing matrix '" << name;
	out << "' (" << this->Rows() << " x " << this->Cols() << ")  Bandwith = "<<fBand<<"\n";
	TPZMatrix<TVar>::Print(0,out,form);
}

/** @brief Overload << operator to output entries of TPZSBMatrixLapack matrix ***/
//ok
template<class TVar>
std::ostream&
operator<<(std::ostream& out,TPZSBMatrixLapack<TVar>  &A)
{
	out.width( 8 );
	out.precision( 4 );
	
	out <<"\n(" << A.Rows() << " x " << A.Cols()
	<< ")  Bandwith = "<< A.GetBand()<<"\n";
	
	for ( long row = 0; row < A.Rows(); row++)
    {
		out << "\t";
		for ( long col = 0; col < A.Cols(); col++ )
			out << A.GetVal( row, col) << "  ";
		out << "\n";
    }
	
	return  out << "\n";
}

/******** Operacoes com matrizes BANDA SIMETRICA  ********/

/******************/
/*** Operator = ***/
//ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator=(const TPZSBMatrixLapack<TVar> &A )
{
	Clear();
	Copy( A );
	return( *this );
}

/******************/
/*** Operator + ***/
//ok
template<class TVar>
TPZSBMatrixLapack<TVar>
TPZSBMatrixLapack<TVar>::operator+(const TPZSBMatrixLapack<TVar> &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrixLapack ) <incompatible dimensions>" );

	TVar *elemMax;
	TVar *elemMin;
	long  sizeMax;
	long  sizeMin;
	if ( fBand > A.fBand )
    {
		sizeMax = fBand;
		elemMax = fDiag;
		sizeMin = A.fBand;
		elemMin = A.fDiag;
    }
	else
    {
		sizeMax = A.fBand;
		elemMax = A.fDiag;
		sizeMin = fBand;
		elemMin = fDiag;
    }
  const long diffSize = sizeMax - sizeMin;
  const long dim = this->Dim();
	TPZSBMatrixLapack<TVar> res( this->Dim(), sizeMax );
	
	TVar *dest = res.fDiag;
	long i;
  for (  i = 0; i < diffSize * dim; i++ )
    *dest++ = *elemMax++;
  for ( ; i < ( sizeMax + 1 ) * dim; i++ )
    *dest++ = (*elemMax++) + (*elemMin++);
	
	return( res );
}

/******************/
/*** Operator - ***/
//ok
template<class TVar>
TPZSBMatrixLapack<TVar>
TPZSBMatrixLapack<TVar>::operator-(const TPZSBMatrixLapack<TVar> &A ) const
{
  if ( this->Dim() != A.Dim() )
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator-( TPZSBMatrixLapack ) <incompatible dimensions>" );
  
  TVar *elemMax;
  TVar *elemMin;
  long  sizeMax;
  long  sizeMin;
  int isThisTheBiggest;
  if ( fBand > A.fBand )
  {
    sizeMax = fBand;
    elemMax = fDiag;
    sizeMin = A.fBand;
    elemMin = A.fDiag;
    isThisTheBiggest = 2;
  }
  else
  {
    sizeMax = A.fBand;
    elemMax = A.fDiag;
    sizeMin = fBand;
    elemMin = fDiag;
    isThisTheBiggest = 0;
  }
  const long diffSize = sizeMax - sizeMin;
  const long dim = this->Dim();
  TPZSBMatrixLapack<TVar> res( this->Dim(), sizeMax );
  
  // SUBTRACT
  TVar *dest = res.fDiag;
  long i;
  for (  i = 0; i < diffSize * dim; i++ )
    *dest++ = (*elemMax++) * real( isThisTheBiggest - 1 );
  for ( ; i < ( sizeMax + 1 ) * dim; i++ )
    *dest++ = (*elemMax++) * real( isThisTheBiggest - 1 ) + (*elemMin++) * real( 1 - isThisTheBiggest );
  
  return( res );
}

/*******************/
/*** Operator += ***/

template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator+=(const TPZSBMatrixLapack<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator+=( TPZSBMatrixLapack ) <incompatible dimensions>" );
	
	// No caso de as bandas serem iguais (ou se "tornarem" iguais).
	if ( fBand <= A.fBand )
    {
		if ( fBand < A.fBand )
			SetBand( A.fBand );
		TVar *pm = fDiag;
		TVar *pa = A.fDiag;
		TVar *end = &fDiag[ Size() ];
		while ( pm < end )
			*pm++ += *pa++;
    }
	else
    {
		// Se a banda desta matriz for maior...
		TVar *pThis  = fDiag;
		TVar *pA     = A.fDiag;
		long sizeThis = fBand + 1;
		long sizeA    = A.fBand + 1;
		long inc = sizeThis - sizeA;
		for ( long col = 0; col < this->Dim(); col++, pThis += inc )
			for ( long i = 0; i < sizeA; i++ )
				*pThis++ += *pA++;
    }
	
	this->fDecomposed = 0;
	return( *this );
}

/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator-=(const TPZSBMatrixLapack<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-=( TPZSBMatrixLapack ) <incompatible dimensions>" );
	
	// No caso de as bandas serem iguais (ou se "tornarem" iguais)
	if ( fBand <= A.fBand )
    {
		if ( fBand < A.fBand )
			SetBand( A.fBand );
		TVar *pm = fDiag;
		TVar *pa = A.fDiag;
		TVar *end = &fDiag[ Size() ];
		while ( pm < end )
			*pm++ -= *pa++;
    }
	else
    {
		// Se a banda desta matriz for maior...
		TVar *pThis  = fDiag;
		TVar *pA     = A.fDiag;
		long sizeThis = fBand + 1;
		long sizeA    = A.fBand + 1;
		long inc = sizeThis - sizeA;
		for ( long col = 0; col < this->Dim(); col++, pThis += inc )
			for ( long i = 0; i < sizeA; i++ )
				*pThis++ -= *pA++;
    }
	
	this->fDecomposed = 0;
	return( *this );
}

template<class TVar>
void TPZSBMatrixLapack<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
								const TVar alpha,const TVar beta ,const int opt,const int stride ) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	if ((!opt && this->Cols()*stride != x.Rows()) || this->Rows()*stride != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSBMatrixLapack::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSBMatrixLapack::MultAdd incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt,stride);
	long rows = this->Rows();
	long xcols = x.Cols();
	long ic, r;
	for (ic = 0; ic < xcols; ic++) {
		long begin, end;
		for ( r = 0; r < rows; r++ ) {
			begin = MAX( r - fBand, 0 );
			end   = MIN( r + fBand + 1, rows );
			TVar val = 0.;
			// Calcula um elemento da resposta.
			for ( long i = begin ; i < end; i++ ) val += GetVal( r, i ) * x.GetVal(stride * i, ic );
			val *= alpha;
			val += z.GetVal(r*stride,ic);
			z.PutVal( r*stride , ic, val );
		}
	}
}

/******** Operacoes com MATRIZES GENERICAS ********/

// Estas operacoes com matrizes genericas, usam a parte triangular
// inferior da maior matriz quadrada de A. Ex.:
//
//  Se A = 01 02 03 04   A matriz usada sera':  01 05 09
//         05 06 07 08                          05 06 10
//         09 10 11 12                          09 10 11
//

/******** Operacoes com valores NUMERICOS ********/
//
// As operacoes com valores numericos sao efetuadas apenas nos
// elementos alocados. Em especial, as operacoes A = 0.0 e A *= 0.0
// desalocam todos os elementos da matriz.
//

/*****************************/
/*** Operator * ( REAL ) ***/

template<class TVar>
TPZSBMatrixLapack<TVar>
TPZSBMatrixLapack<TVar>::operator*(const TVar value ) const
{
	TPZSBMatrixLapack<TVar> res( this->Dim(), fBand );
	
	TVar *pr  = res.fDiag;
	TVar *pm  = fDiag;
	TVar *end = &fDiag[Size()];
	while ( pm < end )
		*pr++ = (*pm++) * value;
	return( res );
}

/******************************/
/*** Operator += ( REAL ) ***/

/******************************/
/*** Operator *= ( REAL ) ***/

template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator*=(const TVar value )
{
	TVar *pm  = fDiag;
	TVar *end = &fDiag[Size()];
	while ( pm < end )
		*pm++ *= value;
	
	this->fDecomposed = 0;
	return( *this );
}

/**************/
/*** Resize ***/
//
// Muda as dimensoes da matriz, mas matem seus valores antigos. Novas
// posicoes sao criadas com ZEROS.
//
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Resize(const long newDim ,const long)
{
	if ( newDim == this->Dim() )
		return( 1 );
	
	if (fBand>this->Dim()-1) fBand=this->Dim()-1;//misael:19/10/96
	// Cria nova matrix.
	long  newSize  = newDim * (fBand + 1);
	TVar *newDiag = new TVar[newSize] ;
	
	// Copia os elementos para a nova matriz.
	TVar *src = fDiag;
	TVar *dst = newDiag;
	TVar *end = &newDiag[ MIN(Size(), newSize) ];
	while ( dst < end )
		*dst++ = *src++;
	
	// Zera as posicoes que sobrarem (se sobrarem).
	end = &newDiag[newSize];
	while ( dst < end )
		*dst++ = 0.0;
	
	// Descarta a matriz antiga e valida a nova matriz.
	if ( fDiag != NULL )
		delete( fDiag );
	fDiag = newDiag;
	this->fCol = this->fRow = newDim;
	this->fDecomposed = 0;
	return( 1 );
}

/*************/
/*** Redim ***/
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Redim(const long newDim ,const long)
{
	if ( newDim != this->Dim() )
    {
		// Aloca a nova matriz.
		this->fRow = this->fCol = newDim;
		if ( this->fDiag != NULL )
			delete( this->fDiag );
		this->fDiag = new TVar[Size()];
		if ( fDiag == NULL )
			TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>" );
    }
	
	TVar *dst = fDiag;
	TVar *end = &fDiag[Size()];
	while ( dst < end )
		*dst++ = 0.;
	
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}

template<class TVar>
int
TPZSBMatrixLapack<TVar>::Zero()
{
	TVar *dst =this->fDiag;
	TVar *end = &fDiag[this->Size()];
	while ( dst < end )
		*dst++ = 0.;
	
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}


/****************/
/*** Set Band ***/
template<class TVar>
int
TPZSBMatrixLapack<TVar>::SetBand(const long newBand )
{
	if ( this->fBand == newBand )
		return( 1 );
	
	if ( this->fBand > (this->Dim() - 1) )
		return( 0 );
	
	TVar *newDiag = new TVar[this->Dim() * (newBand + 1)];
	
	// Copia os elementos antigos para a nova alocacao.
	TVar *pNew = newDiag;
	TVar *pOld = fDiag;
	long newSize  = newBand + 1;
	long oldSize  = fBand + 1;
	long minSize  = MIN( newSize, oldSize );
	long i;
	for ( i = 0; i < minSize; i++ )
    {
		TVar *end = pNew + (this->Dim() * newSize);
		TVar *dst = pNew++;
		TVar *src = pOld++;
		for ( ; dst < end; src += oldSize, dst += newSize )
			*dst = *src;
    }
	
	// Preenche com zero os elementos que sobrarem (se sobrarem).
	for ( ; i < newSize; i++ )
    {
		TVar *end = pNew + (this->Dim() * newSize);
		TVar *dst = pNew++;
		for ( ; dst < end; dst += newSize )
			*dst = 0.0;
    }
	
	if ( fDiag == NULL )
		delete( fDiag );
	fDiag = newDiag;
	
	if ( newBand < fBand )
		this->fDecomposed = 0;
	fBand = newBand;
	
	return( 1 );
}

/********************* Resolucao de sistemas *********************/

/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_Cholesky(std::list<long> &singular)
{
	return Decompose_Cholesky();
}

template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_Cholesky()
{
	if (  this->fDecomposed && this->fDecomposed != ECholesky )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	if (this->fDecomposed) {
		return 1;
	}
	long size = fBand + 1; // Tamanho da banda.
	long next = fBand + 1; // O quanto se soma para "andar" na diagonal.
	TVar *row_k   = fDiag;
	TVar *row_end = fDiag + Size();
	
	for ( long k = 0; k <this-> Dim(); k++, row_k += next  )
    {
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		TVar sum = 0.0;
		TVar *elem_k = row_k + 1;
		TVar *end_k  = row_k + next;
		for ( ; elem_k < end_k; elem_k++ )
			sum += (*elem_k) * (*elem_k);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
      if ( std::abs(*row_k -= sum) < 1.e-10 )
			return( 0 );
		*row_k = sqrt( *row_k );
		
		// Loop para i = k+1 ... Dim().
		//
		TVar *row_i = row_k + next;
		for ( long j = 2; row_i < row_end; j++, row_i += next )
		{
			// Se tiverem elementos na linha 'i' cuja coluna e'
			//  menor do que 'K'...
			if ( j < size )
			{
				// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
				sum = 0.0;
				TVar *elem_i = row_i + j;
				TVar *end_i  = row_i + size;
				elem_k = row_k + 1;
				while ( elem_i < end_i )
					sum += (*elem_i++) * (*elem_k++);
				
				// Faz A(i,k) = (A(i,k) - sum) / A(k,k)
				row_i[j-1] = (row_i[j-1] -sum) / (*row_k);
			}
			
			// Se nao tiverem estes elementos, sum = 0.0.
			else if ( j == size )
				row_i[j-1] /= (*row_k);
			
			// Se nao existir nem o elemento A(i,k), nao faz nada.
		}
    }
	
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}

/**********************/
/*** Decompose LDLt ***/
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_LDLt(std::list<long> &singular)
{
	return Decompose_LDLt();
}

template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_LDLt()
{
	
	if (  this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
	long j,k,l, begin,end;
	TVar sum;
	
	
	for ( j = 0; j < this->Dim(); j++ )
    {
		//Print("curernt");
		sum=0.;
		
		begin = MAX( long(j - fBand), 0 );
		//cout<<"begin="<<begin<<"\n";
		for ( k=begin; k<j; k++)
		{
			sum=sum-GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
			//cout<<"(k,j)"<<k<<" "<<j<<"\n";
		}
		
		
		//	 operator()(j,j)=GetVal(j,j)+sum;
		PutVal(j,j,GetVal(j,j)+sum);
		//cout<<"\n(j,j)"<<j<<" "<<j<<"\n\n";
		for ( k=0; k<j; k++)
		{
			end   = MIN( long(k + fBand )+1, this->Dim() );
			for( l=j+1; l<end;l++)
			{
				PutVal(l,j, GetVal(l,j)-GetVal(k,k)*GetVal(j,k)*GetVal(l,k) );
				/*cout<<"end="<<end<<"\n";
				 cout<<"(l,j)"<<l<<" "<<j<<"\n";
				 cout<<"(j,k)"<<j<<" "<<k<<"\n";
				 cout<<"(l,k)"<<l<<" "<<k<<"\n\n";
				 */
			}
		}
		
		if ( IsZero(GetVal(j,j)) ) TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Zero on diagonal>" );
		end  = MIN( long(j + fBand )+1, this->Dim() );
		//cout<<"end="<<end<<"\n";
		for( l=j+1; l<end;l++)
		{
			//cout<<"(l,j)"<<l<<" "<<j<<"\n";
			PutVal( l,j,GetVal(l,j)/GetVal(j,j) ) ;
		}
    }
	this->fDecomposed  = 1;
	this->fDefPositive = 0;
	
	return( 1 );
	
}

/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular inferior.
//
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Subst_Forward( TPZFMatrix<TVar>*B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	
	long size = fBand + 1;
	TVar *row_k = fDiag;
	for ( long k = 0; k < this->Dim(); k++, row_k += size )
		for ( long j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			long end=(k-fBand>0)? fBand:k;//misael
			TVar sum = 0.0;
			TVar *elem_ki = row_k + 1;
			TVar *end_ki  = row_k + end;  //misael
			for ( long i = k-1; elem_ki < end_ki ; i-- )
				sum += (*elem_ki++) * B->GetVal( i, j );
			
			// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, (B->GetVal(k, j) - sum) / *row_k );
		}
	
	return( 1 );
}

template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	long k,j,i,jmax,stepcol=fBand+2;
	for(k=0; k<B->Cols() ; k++)
    {
		for(i=this->Rows()-1; i>=0; i--)
		{
			TVar *diagk=fDiag+(fBand+1)*i;
			jmax=( (i+fBand)>this->Rows()-1)? this->Rows()-1 : i+fBand;
			for(j=i+1;j<=jmax;j++)
			{
				diagk+=stepcol;
				B->operator()(i,k)-=B->GetVal(j,k)*(*diagk);
			}
			B->operator()(i,k)/=GetVal(i,i);
			
		}
		
    }
	
	return ( 1 ) ;
	
}

/***********************/
/*** Subst L Forward ***/
//
//  Faz a "Forward substitution" assumindo que os elementos
//   da diagonal sao todos iguais a 1.
//
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LForward-> uncompatible matrices") ;
	
	long size = fBand + 1;
	TVar *row_k = fDiag;
	long i,j,k;
	for ( k = 0; k < this->Dim(); k++, row_k += size )
		for ( j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			
			long end=(k-fBand>0)? fBand:k;  //misael
			TVar sum = 0.0;
			TVar *elem_ki = row_k + 1;
			TVar *end_ki  = row_k + end;
			for ( i = k-1; elem_ki <= end_ki ; i-- )//misael
				sum += (*elem_ki++) * B->GetVal( i, j );
			
			// Faz b[k] = (b[k] - sum).
			//
			B->PutVal( k, j, (B->GetVal( k, j ) - sum) );
			
		}
	
	return( 1 );
}

/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const
{
	
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Diag-> uncompatible matrices") ;
	
	
	long size = fBand + 1;
	TVar *row_k = fDiag;
	for ( long k = 0; k < this->Dim(); k++, row_k += size )
		for ( long j = 0; j < B->Cols(); j++ )
			B->PutVal( k, j, B->GetVal( k, j) / *row_k );
	
	return( 1 );
}

template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LBackward-> uncompatible matrices") ;
	
	long k,j,i,jmax,stepcol=fBand+2;
	
	for(k=0; k<B->Cols() ; k++)
    {
		for(i=this->Rows()-1; i>=0; i--)
		{
			TVar *diagk=fDiag+(fBand+1)*i;
			jmax=( (i+fBand)>this->Rows()-1)? this->Rows()-1 : i+fBand;
			for(j=i+1;j<=jmax;j++)
			{
				diagk+=stepcol;
				B->operator()(i,k)-=B->GetVal(j,k)*(*diagk);
			}
		}
		
    }
	
	return 1;
	
}

/**************************** PRIVATE ****************************/

/*************/
/*** CLear ***/
//ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Clear()
{
	if ( this->fDiag != NULL )
		delete []fDiag ;
	this->fRow = this->fCol = 0;
	fDiag = NULL;
	this->fDecomposed = 0;
	return( 1 );
}

/************/
/*** Copy ***/
//ok
template<class TVar>
void
TPZSBMatrixLapack<TVar>::Copy(const TPZSBMatrixLapack<TVar> &A )
{
	this->fBand = A.fBand;
	this->fRow = this->fCol = A.Dim();
	this->fDiag = new TVar[Size()];
	this->fDecomposed  = A.fDecomposed;
	this->fDefPositive = A.fDefPositive;
	
	if ( fDiag == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Copy( TPZSBMatrixLapack ) <memory allocation error>" );
	
	TVar *dst = fDiag;
	TVar *src = A.fDiag;
	TVar *end = A.fDiag + Size();
	while ( src < end )
		*dst++ = *src++;
}

#ifdef OOPARLIB

template<class TVar>
int TPZSBMatrixLapack<TVar>::Unpack (TReceiveStorage *buf ){
	TSimMatrix::Unpack(buf);
	buf->UpkInt(&fBand);
	buf->UpkDouble(fDiag,fBand);
	return 1;
}


template<class TVar>
TSaveable *TPZSBMatrixLapack<TVar>::Restore(TReceiveStorage *buf) {
	TPZSBMatrixLapack<TVar> *m = new TPZSBMatrixLapack<TVar>();
	m->Unpack(buf);
	return m;
}

template<class TVar>
int TPZSBMatrixLapack<TVar>::Pack( TSendStorage *buf ){
	TSimMatrix::Pack(buf);
	buf->PkInt(&fBand);
	buf->PkDouble(fDiag,fBand);
	return 1;
}

template<class TVar>
int TPZSBMatrixLapack<TVar>::DerivedFrom(long Classid){
	if(Classid == GetClassID()) return 1;
	return TSimMatrix::DerivedFrom(Classid);
}

template<class TVar>
int TPZSBMatrixLapack<TVar>::DerivedFrom(char *classname){
	if(!strcmp(ClassName(),classname)) return 1;
	return TSimMatrix::DerivedFrom(classname);
}

#endif

// Inicializando os templates
template class TPZSBMatrixLapack< double >;
template class TPZSBMatrixLapack< long double >;
template class TPZSBMatrixLapack< complex<double> >;

