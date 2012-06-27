/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "pzcmesh.h"
#include "pzelchdiv.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "pzmaterialdata.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.TPZCompElHDiv"));
#endif

using namespace std;


template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, int &index) :
TPZIntelGen<TSHAPE>(mesh,gel,index,1) {
	int i;
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	int nconflux= TPZCompElHDiv::NConnects();
    this->fConnectIndexes.Resize(nconflux);
	gel->SetReference(this);
	
	for(i=0;i< nconflux;i++)
	{
        int sideaux= i + TSHAPE::NCornerNodes;
		this->fConnectIndexes[i] = this->CreateMidSideConnect(sideaux);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "After creating last flux connect " << i << std::endl;
            //	this->Print(sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
        
		mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
		this->IdentifySideOrder(sideaux);
    }	
    
	
    int sideorder = SideOrder(TSHAPE::NSides-1);
	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	TPZManVector<int,3> order(3,sideorder);
	this->fIntRule.SetOrder(order);
}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh, const TPZCompElHDiv<TSHAPE> &copy) :
TPZIntelGen<TSHAPE>(mesh,copy)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		this-> fConnectIndexes[i] = copy.fConnectIndexes[i];
	}
}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv(TPZCompMesh &mesh,
									 const TPZCompElHDiv<TSHAPE> &copy,
									 std::map<int,int> & gl2lcConMap,
									 std::map<int,int> & gl2lcElMap) :
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		int lcIdx = -1;
		int glIdx = copy.fConnectIndexes[i];
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			this-> fConnectIndexes[i] = -1;
			return;
		}
		this-> fConnectIndexes[i] = lcIdx;
	}
	
}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::TPZCompElHDiv() :
TPZIntelGen<TSHAPE>()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		this-> fConnectIndexes[i] = -1;
	}
}

template<class TSHAPE>
TPZCompElHDiv<TSHAPE>::~TPZCompElHDiv(){
	
}

template<class TSHAPE>
MElementType TPZCompElHDiv<TSHAPE>::Type() {
	return TSHAPE::Type();
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NConnects() const {
	int dimension = Dimension()-1;
	
	return TSHAPE::NumSides(dimension) + 1;
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetConnectIndex(int i, int connectindex){
#ifndef NODEBUG
	if(i<0 || i>= this->NConnects()) {
		std::cout << " TPZCompElHDiv<TSHAPE>::SetConnectIndex index " << i <<
		" out of range\n";
		DebugStop();
		return;
	}
#endif
	this-> fConnectIndexes[i] = connectindex;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NConnectShapeF(int connect)const
{
    if(connect >= this->NConnects())
    {
        PZError << "TPZCompElHDiv<TSHAPE>::NConnectShapeF: there is not this connect " <<  endl;
        return -1;
    }
    
    if (TSHAPE::Type()==EQuadrilateral)
    {
        int iside = connect+TSHAPE::NCornerNodes;
        if(TSHAPE::SideDimension(iside)< this->Dimension()-2)
        {
            PZError << "TPZCompElHDiv<TSHAPE>::NConnectShapeF: no shape associate " <<  endl;
            return -1;
                
        }
        int order = ConnectOrder(connect);
        if(order < 0) return 0;
        
        TPZStack<int> smallsides;
        TSHAPE::LowerDimensionSides(iside,smallsides);
        //i.e., trata-se do lado de mesma dimensao que o elemento
        if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension)
        {
            int NShapeFace = 0;
            //funcoes para os lados menor que o proprio elemento
            for(int nside=TSHAPE::NCornerNodes; nside<smallsides.NElements();nside++)
            {
                int sideorder= this->SideOrder(nside);
                NShapeFace += TSHAPE::NConnectShapeF(nside,sideorder);
            }
            //se a ordem for maior q um tiraremos a ultima função interna
            if (order>1)
            {
                int aux=TSHAPE::NConnectShapeF(iside,order)-1;//-1;
                int nshape=NShapeFace+ TSHAPE::Dimension*(aux)-4*(order-2);
#ifdef LOG4CXX
                {
                    std::stringstream sout;
                    sout << endl<<" Connect : " << connect << " numero de func de forma " << nshape<<std::endl;
                    LOGPZ_DEBUG(logger, sout.str().c_str());
                }
#endif
                return(nshape);//pois nem todas as funcoes levam as duas direcoes 
            }
            else {
                int aux=TSHAPE::NConnectShapeF(iside,order);
                return(NShapeFace+ TSHAPE::Dimension*(aux));
            }
                
        }
            
        //trata-se do lado de 1 dimensao a menos que a dimensao do elemento
        else if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension-1) 
        {
            int NShapeF = 0;
            for(int j=0;j< smallsides.NElements();j++)
            {
                NShapeF += TSHAPE::NConnectShapeF(j,order-1);
            }
                
            //se a ordem for maior q um tbem tira-se a ultima funcao de cada lado
            if (order>1) {
                int result=NShapeF + TSHAPE::NConnectShapeF(iside,order)-1;//era -1..coloquei -2
#ifdef LOG4CXX
            {
                std::stringstream sout;
                sout << endl<<" Connect : " << connect << " numero de func de forma " << result<<std::endl;
                LOGPZ_DEBUG(logger, sout.str().c_str());
            }
#endif
                return(result);
                    
            }
            else {
                int result=NShapeF + TSHAPE::NConnectShapeF(iside,order);
                return(result);
            }
                
        }
    }	
    
    if (TSHAPE::Type()==ETriangle){
        
        if(connect< NConnects()){//tirando o connect da pressao
            //int iside = TSHAPE::NSides - NConnects() + connect+1;//é o lado do elemento q estou, acrescentei um connect a mais
            int iside = connect+TSHAPE::NCornerNodes;
            if(TSHAPE::SideDimension(iside)< this->Dimension()-2)
            {
                PZError << "TPZCompElHDiv<TSHAPE>::NConnectShapeF: no shape associate " <<  endl;
                return -1;
                
            }
            int order = ConnectOrder(connect);
            
            if(order < 0) return 0;
            
            TPZStack<int> smallsides;
            TSHAPE::LowerDimensionSides(iside,smallsides);
            if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension)//i.e., trata-se do lado de mesma dimensao que o elemento
            {
                int NShapeFace = 0;
                for(int nside=TSHAPE::NCornerNodes; nside<smallsides.NElements();nside++){
                    int sideorder= this->SideOrder(nside);
                    NShapeFace += TSHAPE::NConnectShapeF(nside,sideorder);
                    
                }
                return (NShapeFace + TSHAPE::Dimension*TSHAPE::NConnectShapeF(iside,order))  ;
                
            }
            else if(TSHAPE::SideDimension(iside) == TSHAPE::Dimension-1)//i.e., trata-se do lado de 1 dimensao a menos que a dimensao do elemento
            {
                int NShapeF = 0;
                for(int j=0;j< smallsides.NElements();j++)
                {
                    NShapeF += TSHAPE::NConnectShapeF(j,order);
                }
                int result=NShapeF + TSHAPE::NConnectShapeF(iside,order);
                return (result);
            }
        }
    }
    
    else {
        std::cout<< "No implemented yet"<<std::endl;
        DebugStop();
    }
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout <<__PRETTY_FUNCTION__<< "unhandled case ";
        LOGPZ_ERROR(logger,sout.str())
    }
#endif
    return -1;
    
    
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetIntegrationRule(int ord) {
	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	this->fIntRule.SetOrder(order);
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NSideConnects(int side) const{
	if(TSHAPE::SideDimension(side)<= Dimension()-2) return 0;
	if(TSHAPE::SideDimension(side)==Dimension()-1) return 1;
	if(TSHAPE::SideDimension(side)== Dimension()) return this->NConnects();
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case ";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	return -1;
	
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::SideConnectLocId(int node,int side) const {
	if(TSHAPE::SideDimension(side)<= TSHAPE::Dimension - 2 || node >= NSideConnects(side)) {
		PZError << "TPZCompElHDiv<TSHAPE>::SideConnectLocId no connect associate " <<  endl;
		return -1;
	}
	
    return side-TSHAPE::NCornerNodes; 
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ConnectSideLocId(int connect) const{
    
    int side = connect + TSHAPE::NCornerNodes ; 
    return side;   
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(NConnects());
	int i;
	for(i=0; i<NConnects(); i++) {
		ord[i] = ConnectOrder(i);
	}
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::PreferredSideOrder(int side) {
	if(TSHAPE::SideDimension(side) < Dimension()-1)
	{
		PZError << __PRETTY_FUNCTION__ << " side " << side << std::endl;
	}
	int connect= SideConnectLocId(0,side);
	if(connect<0 || connect > NConnects()) {
		PZError << "TPZCompElHDiv<TSHAPE>::PreferredSideOrder no polynomial associate " <<  endl;
		return -1;
	}
	if(connect<NConnects()) {
		int order =this->fPreferredOrder;
		return this->AdjustPreferredSideOrder(side,order);
	}
	PZError << "TPZCompElHDiv<TSHAPE>::PreferredSideOrder called for connect = " << connect << "\n";
	return 0;
	
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ConnectIndex(int con) const{
#ifndef NODEBUG
	if(con<0 || con>= this->NConnects()) {
		std::cout << "TPZCompElHDiv::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << this-> NConnects() << std::endl;
		DebugStop();
		return -1;
	}
	
#endif
	
	return this->fConnectIndexes[con];
}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetPreferredOrder(int order)
{
	this->fPreferredOrder = order;
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SetSideOrder(int side, int order){
	int connectaux= SideConnectLocId(0,side);
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElHDiv::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	}
	TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order);
    int seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial * mat =this-> Material();
    if(mat) nvar = mat->NStateVariables();
    c.SetNState(nvar);
    int nshape = NConnectShapeF(connectaux);
    c.SetNShape(nshape);
	this-> Mesh()->Block().Set(seqnum,nshape*nvar);
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::ConnectOrder(int connect) const{
	if (connect < 0 || connect >= this->NConnects()){
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Connect index out of range connect " << connect <<
			" nconnects " << NConnects();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return -1;
	}
	
	if (this->fConnectIndexes[connect] == -1) {
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " connect " << connect
		<< " is not initialized" << std::endl;
#ifdef LOG4CXX
		LOGPZ_ERROR(logger,sout.str());
#else
		std::cout << sout.str() << std::endl;
#endif
		return 0;
	}
	
    TPZConnect &c = this-> Connect(connect);
    return c.Order();
}

template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::SideOrder(int side) const
{
	if(!NSideConnects(side)) return -1;
	int corder =SideConnectLocId(0, side);
	int maxorder = 0;
	int conectaux;
	if(corder>=0 || corder <= NConnects()) return ConnectOrder(corder);
    
	TPZStack< int > high;
	TSHAPE::HigherDimensionSides(side, high);
	int highside= high.NElements();
	
	
	for(int j=0;j<highside;j++)
	{
		conectaux =SideConnectLocId(0, high[j]);
		maxorder = (ConnectOrder(conectaux) > maxorder) ? ConnectOrder(conectaux) : maxorder;
	}
	
	return maxorder;
	
}

/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::FirstShapeIndex(TPZVec<int> &Index){
	
	Index.Resize(TSHAPE::NSides+1);
	Index[0]=0;
	
	for(int iside=0;iside<TSHAPE::NSides;iside++)
	{
		int order= SideOrder(iside);
		Index[iside+1] = Index[iside] + TSHAPE::NConnectShapeF(iside,order);
		
	}
	
#ifdef LOG4CXX
    std::stringstream sout;
    sout << "First  Index" << Index;
    LOGPZ_DEBUG(logger,sout.str())
#endif
}


template<class TSHAPE>
int TPZCompElHDiv<TSHAPE>::NFluxShapeF() const{
    int in,result=0;
    int nn=TPZCompElHDiv::NConnects();
    for(in=0;in<nn;in++) result += NConnectShapeF(in);
    return result;
}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int> > & ShapeAndVec){
    
    
    // VectorSide indicates the side associated with each vector entry
    TPZManVector<int,27> FirstIndex;
    // the first index of the shape functions
    FirstShapeIndex(FirstIndex);
    
    int count=0;
    int nshapeflux= NFluxShapeF();
    ShapeAndVec.Resize(nshapeflux);
    
    if (TSHAPE::Type()==EQuadrilateral) {
        
        TPZManVector<int,4> ids(4,0);
        TPZGeoEl *gel = this->Reference();
        int id;
        for (id=0; id<4; id++) {
            ids[id] = gel->NodePtr(id)->Id();
        }
        
        for(int jvec=0;jvec< VectorSide.NElements();jvec++)
            
        {
            if (jvec==2||jvec==5||jvec==8||jvec==11){
                int lside=VectorSide[jvec];
                int fshape1= FirstIndex[lside];
                int fshape2= FirstIndex[lside+1]-1;
                for (int ishape=fshape1; ishape<fshape2; ishape++)//estou tentando tirar a ultima funcao
                {
#ifdef LOG4CXX
                    std::stringstream sout;
                    sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                    LOGPZ_DEBUG(logger,sout.str())
#endif
                    ShapeAndVec[count++]=std::pair<int,int>(jvec,ishape);
                    
                }
                
                
            }
            else if(jvec==16 || jvec ==17)
            {
                int lside = VectorSide[jvec];
                int order = SideOrder(lside);
                int nshape = TSHAPE::NConnectShapeF(lside,order-1);//order-1);
                TPZFNMatrix<25> sideorders(2,nshape);
                int ksi,eta;
                int loccount = 0;
                int transid = TSHAPE::GetTransformId(ids);
                switch (transid) 
                {
                    case 0:
                    case 3:
                    case 4:
                    case 7:
                        for (ksi=0; ksi<order-2; ksi++) {
                            for (eta=0; eta<order-2; eta++) {
                                sideorders(0,loccount) = ksi+2;
                                sideorders(1,loccount) = eta+2;
                                loccount++;
                            }
                        }
                        break;
                    case 2:
                    case 6:
                    case 1:
                    case 5:
                        for (eta=0; eta<order-2; eta++) {
                            for (ksi=0; ksi<order-2; ksi++) {
                                sideorders(0,loccount) = ksi+2;
                                sideorders(1,loccount) = eta+2;
                                loccount++;
                            }
                        }
                        
                        break;
                    default:
                        DebugStop();
                        break;
                }
                int ish=0;
                int fshape1 = FirstIndex[lside];
#ifdef LOG4CXX
                std::stringstream sout;
                sideorders.Print("SideOrders= ", sout ,EFormatted);
                LOGPZ_DEBUG(logger,sout.str())
#endif
                for (ish=0; ish<nshape; ish++) {//coloquei -2
                    if (jvec ==16 && sideorders(0,ish) <= order && sideorders(1,ish) < order) {
#ifdef LOG4CXX
                        std::stringstream sout;
                        sout << " <vec,shape> " << "< "<<jvec << " * "<<fshape1+ish << "> "<<std::endl;
                        sout << " side order ksi " << sideorders(0,ish) << " side order eta " << sideorders(1,ish);
                        LOGPZ_DEBUG(logger,sout.str())
#endif
                        ShapeAndVec[count++]=std::pair<int,int>(jvec,fshape1+ish);
                    }
                    if (jvec ==17 && sideorders(0,ish) < order && sideorders(1,ish) <= order) {
#ifdef LOG4CXX
                        std::stringstream sout;
                        sout << " <vec,shape> " << "< "<<jvec << " * "<<fshape1+ish << "> "<<std::endl;
                        sout << " side order ksi " << sideorders(0,ish) << " side order eta " << sideorders(1,ish);
                        LOGPZ_DEBUG(logger,sout.str())
#endif
                        ShapeAndVec[count++]=std::pair<int,int>(jvec,fshape1+ish);
                    }
                }
            }
            
            else{
                
                
                int lside=VectorSide[jvec];
                int fshape1= FirstIndex[lside];
                int fshape2= FirstIndex[lside+1];
                for (int ishape=fshape1; ishape<fshape2; ishape++){
#ifdef LOG4CXX
                    std::stringstream sout;
                    sout << " <vec,shape> " << "< "<<jvec << " * "<<ishape << "> "<<std::endl;
                    LOGPZ_DEBUG(logger,sout.str())
#endif
                    ShapeAndVec[count++]=std::pair<int,int>(jvec,ishape);
                    
                }
            }
            
        }
        
    }
    
    else{
        
        int count=0;
        for(int jvec=0;jvec< VectorSide.NElements();jvec++)
        {
            int lside=VectorSide[jvec];
            int fshape1= FirstIndex[lside];
            int fshape2= FirstIndex[lside+1];
            for (int ishape=fshape1; ishape<fshape2; ishape++)
            {
                ShapeAndVec[count++]=std::pair<int,int>(jvec,ishape);
            }
            
        }
        
    }
    
    
#ifdef LOG4CXX
    std::stringstream sout;
    sout << " ShapeAndVec" << ShapeAndVec;
    LOGPZ_DEBUG(logger,sout.str())
#endif
    
}


//compute the values of the shape function of the side
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
    if(side==TSHAPE::NSides){
        std::cout<<"Don't have side shape associated to this side";
        DebugStop();
    }
	if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension -1 ){
		return ;
	}
    int ncontained = TSHAPE::NContainedSides(side);
    int nsideshape = 0;
    int order = SideOrder(side);
    int is;
    for (is=0; is<ncontained; is++) {
        int ic = TSHAPE::ContainedSideLocId(side,is);
        nsideshape += TSHAPE::NConnectShapeF(ic,order);
    }
    if (nsideshape != this->NSideShapeF(side)) {
        // create a philoc and copy
        TPZFNMatrix<200> philoc(nsideshape,1),dphiloc(TSHAPE::SideDimension(side),nsideshape);
        TPZIntelGen<TSHAPE>::SideShapeFunction(side,point,philoc,dphiloc);
        int nsh = phi.Rows();
        for (int ish = 0; ish < nsh; ish++) {
            phi(ish,0) = philoc(ish,0);
            for (int d=0; d<TSHAPE::SideDimension(side); d++) {
                dphi(d,ish) = dphiloc(d,ish);
            }
        }
    }
    else
    {
        TPZIntelGen<TSHAPE>::SideShapeFunction(side,point,phi,dphi);
    }
    
    if(TSHAPE::SideDimension(side) == 1)
    {
        int locnod0 = TSHAPE::SideNodeLocId(side,0);
        int locnod1 = TSHAPE::SideNodeLocId(side,1);
        TPZGeoEl *gel = this->Reference();
        int locnodid0 = gel->NodePtr(locnod0)->Id();
        int locnodid1 = gel->NodePtr(locnod1)->Id();
        if(locnodid0 > locnodid1)
        {
            REAL temp;
            temp = phi(0,0);
            phi(0,0) = phi(1,0);
            phi(1,0) = temp;
            int d;
            for (d=0; d<TSHAPE::SideDimension(side); d++) {
                temp = dphi(d,0);
                dphi(d,0) = dphi(d,1);
                dphi(d,1) = temp;
            }
        }
    }
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol){
	TPZMaterialData data;
	InitMaterialData(data);
    
//	this->ComputeShape(qsi, data.x,data.jacobian,data.axes, data.detjac,data.jacinv,data.phi, data.dphix);
	
	
	this->ComputeSolutionHDiv(data);
	this->Material()->Solution(data,var,sol);
    
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeSolutionHDiv(TPZVec<REAL> &qsi, TPZMaterialData &data){
//	this->ComputeShape(qsi, data.x,data.jacobian,data.axes, data.detjac,data.jacinv,data.phi, data.dphix);
    this->ComputeShape(qsi,data);
    this->ComputeSolutionHDiv(data);
	
    
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                                            const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol){
    TPZMaterialData data;
    data.phi = phi;
    data.dphix = dphix;
    data.axes = axes;
    InitMaterialData(data);
    this->ComputeSolutionHDiv(data);
    sol = data.sol;
    dsol = data.dsol;
	
    
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data){
    
    //this->InitMaterialData(data);
    //this->ComputeShape(qsi,data);
    this->ComputeSolution(qsi, data.phi, data.dphix, data.axes, data.sol, data.dsol);
    
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes) {
	
	TPZGeoEl * ref = this->Reference();
	const int nshape = this->NShapeF();
	const int dim = ref->Dimension();
    
    TPZMaterialData data;
    data.phi.Resize(nshape, 1);
    data.dphix.Resize(dim,nshape);
	data.jacobian.Resize(dim, dim);
    data.jacinv.Resize(dim,dim);
	data.x.Resize(3,0.);
    
	//this->ComputeShape(qsi,x,jacobian,axes,detjac,jacinv,phix,dphix);
    //this->ComputeSolution(qsi, phix, dphix, axes, sol, dsol);
    this->ComputeShape(qsi,data);
    this->ComputeSolution(qsi, data.phi, data.dphix, data.axes, data.sol, data.dsol);
	
	
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeSolutionHDiv(TPZMaterialData &data){
    const int numdof = this->Material()->NStateVariables();
    const int ncon = this->NConnects();
    
	TPZBlock<STATE> &block =this->Mesh()->Block();
    TPZFMatrix<STATE> &MeshSol = this->Mesh()->Solution();
    int numbersol = MeshSol.Cols();
    
	int nsol= this->Dimension()+1;
    data.sol.Resize(numbersol);
    data.dsol.Resize(numbersol);
    for (int is=0; is<numbersol; is++) {
        data.sol[is].Resize(nsol,1);//2 components to the flow  and 1 to the div
        data.sol[is].Fill(0);
		
    }
	//solution associated with the flow
	int iv = 0,ishape=0,ivec=0,cols, jv=0;
    for(int in=0; in<ncon ; in++) {
		TPZConnect *df = &this->Connect(in);
		int dfseq = df->SequenceNumber();
		int dfvar = block.Size(dfseq);
		int pos = block.Position(dfseq);
		
		for(int jn=0; jn<dfvar; jn++) {
			ivec=data.fVecShapeIndex[jv ].first;
			ishape=data.fVecShapeIndex[jv].second;
			
			TPZFNMatrix<3> ivecDiv(3,1);
			ivecDiv(0,0) = data.fNormalVec(0,ivec);
			ivecDiv(1,0) = data.fNormalVec(1,ivec);
			ivecDiv(2,0) = data.fNormalVec(2,ivec);
			TPZFNMatrix<3> axesvec(3,1);
			data.axes.Multiply(ivecDiv,axesvec);
			
			for (int ilinha=0; ilinha<this->Dimension(); ilinha++) {
				cols=iv%numdof;
				for (int is=0; is<numbersol; is++) {
                    data.sol[is][ilinha] += (STATE)data.fNormalVec(ilinha,ivec)*(STATE)data.phi(ishape,0)*MeshSol(pos+jn,is);
                    data.sol[is][nsol-1] += (STATE)axesvec(ilinha,0)*(STATE)data.dphix(ilinha,ishape)*MeshSol(pos+jn,is);//divergente
                    
                }
				
			}
			
			jv++;
		}
		
		iv++;
	}
	
}


template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12)
{
	bool Is_u1PHI = (u1.Cols() == 1) ? true : false;
	bool Is_u2PHI = (u2.Cols() == 1) ? true : false;
	
	if(Is_u1PHI && Is_u2PHI)
	{
		int nu1 = u1.Rows(),nu2 = u2.Rows();
		u12.Redim(nu1+nu2,1);
		int i;
		for(i=0; i<nu1; i++) u12(i,0) = u1(i,0);
		for(i=0; i<nu2; i++) u12(i+nu1,0) = u2(i,0);
		
		
	}
	else if(!Is_u1PHI || !Is_u2PHI) 
	{
		int ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
		int ru12 = ru1 < ru2 ? ru2 : ru1;
		int cu12 = cu1+cu2;
		u12.Redim(ru12,cu12);
		int i,j;
		for(i=0; i<ru1; i++) for(j=0; j<cu1; j++) u12(i,j) = u1(i,j);
		for(i=0; i<ru2; i++) for(j=0; j<cu2; j++) u12(i,j+cu1) = u2(i,j);
	}
	else
	{
		PZError << "TPZCompElHDiv::Append. Bad input parameters " << std::endl;		
		
	}
	
}

template<class TSHAPE>

void TPZCompElHDiv<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	TPZManVector<int,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
	TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
	int i;
	TPZGeoEl *ref = this->Reference();
	for(i=0; i<TSHAPE::NCornerNodes; i++) {
		id[i] = ref->NodePtr(i)->Id();
	}
	
    int nconflux=TPZCompElHDiv::NConnects();
    if(TSHAPE::Type()==EQuadrilateral)
    {
        for(i=0; i< nconflux-1; i++) {
            ord[i] = ConnectOrder(i);
		
            //		 #ifdef LOG4CXX
            //		 std::stringstream sout;
            //		 sout << " ordem do connect " << i << " " << ord[i]<<endl;
            //		 LOGPZ_DEBUG(logger,sout.str())
            //		 #endif
        
        }
        ord[i] = ConnectOrder(i)-1;
    }
    else{
        for(i=0; i< nconflux; i++)
        {
            ord[i] = ConnectOrder(i);
		
        //		 #ifdef LOG4CXX
        //		 std::stringstream sout;
        //		 sout << " ordem do connect " << i << " " << ord[i]<<endl;
        //		 LOGPZ_DEBUG(logger,sout.str())
        //		 #endif
        
        }
    }
	int dimension= TSHAPE::Dimension;
	
	
	const int nshapecont = TSHAPE::NShapeF(ord);
    
    phi.Resize(nshapecont, 1);
    dphi.Resize(dimension, nshapecont);
	TSHAPE::Shape(pt,id,ord,phi,dphi);
    
}


template<class TSHAPE>
TPZTransform TPZCompElHDiv<TSHAPE>::TransformSideToElement(int side){
	return TSHAPE::TransformSideToElement(side);
}

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int> &shapeindex){
	
	TPZManVector<int> firstshapeindex;
	FirstShapeIndex(firstshapeindex);
	int nshape = TPZIntelGen<TSHAPE>::NShapeF();
	shapeindex.Resize(nshape);
	int nsides = sides.NElements();
	int is, count=0;
	for(is=0 ; is<nsides; is++)
	{
		int side = sides[is];
		int sideorder= this->SideOrder(side);
		int NShapeFace = TSHAPE::NConnectShapeF(side,sideorder);
		int ishapeface;
		for(ishapeface=0; ishapeface<NShapeFace; ishapeface++)
		{
			shapeindex[count++] = is;
		}
	}
	shapeindex.Resize(count);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "count = " << count << " nshape " << nshape;
		sout << endl<<"sides associated with the normals "<< sides <<
		"\nnormal associated with each shape function : shape function indexes " << shapeindex;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
	
#ifdef LOG4CXX
	{
		LOGPZ_DEBUG(logger,"Initializing normal vectors")
	}
#endif
	TPZManVector<int> normalsides;
	TPZIntelGen<TSHAPE>::Reference()->ComputeNormals(data.fNormalVec, normalsides);
	// vecindex : lado associado a cada normal
	// vecindex indica apenas o numero do lado associado a cada normal
	// agora temos que expandir para formar pares : vecIndex e shapeindex
	//	ComputeShapeIndex(data.fVecIndex,data.fShapeIndex);
	//data.numberdualfunctions = NConnectShapeF(NConnects()-1);
	IndexShapeToVec(normalsides,data.fVecShapeIndex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		data.fNormalVec.Print("Normal vectors", sout);
		sout << "NormalVector/Shape indexes " << data.fVecShapeIndex << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

// Save the element data to a stream
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
	TPZInterpolatedElement::Write(buf,withclassid);
	TPZManVector<int,3> order(3,0);
	this-> fIntRule.GetOrder(order);
	this->  WriteObjects(buf,order);
	buf.Write(this->fConnectIndexes,TSHAPE::NSides);
	buf.Write(&this->fPreferredOrder,1);
	int classid = this->ClassId();
	buf.Write ( &classid, 1 );
}


// Read the element data from a stream
template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZInterpolatedElement::Read(buf,context);
	TPZManVector<int,3> order;
	this-> ReadObjects(buf,order);
	this-> fIntRule.SetOrder(order);
	buf.Read(this->fConnectIndexes,TSHAPE::NSides);
	buf.Read(&this->fPreferredOrder,1);
	int classid = -1;
	buf.Read( &classid, 1 );
	if ( classid != this->ClassId() )
	{
		std::stringstream sout;
		sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " for an package of id = " << classid;
		LOGPZ_ERROR ( logger, sout.str().c_str() );
	}
}


#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"
#include "tpzpyramid.h"

#include "pzmeshid.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

template<>
void TPZCompElHDiv<TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}



template<class TSHAPE>
void TPZCompElHDiv<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}

template<>
int TPZCompElHDiv<TPZShapePoint>::ClassId() const
{
	return TPZHDIVPOINTID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapePoint>, TPZHDIVPOINTID>;


template<>
int TPZCompElHDiv<TPZShapeLinear>::ClassId() const
{
	return TPZHDIVLINEARID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapeLinear>, TPZHDIVLINEARID>;

template<>
int TPZCompElHDiv<TPZShapeTriang>::ClassId() const
{
	return TPZHDIVTRIANGLEID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapeTriang>, TPZHDIVTRIANGLEID>;

template<>
int TPZCompElHDiv<TPZShapeQuad>::ClassId() const
{
	return TPZHDIVQUADID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapeQuad>, TPZHDIVQUADID>;

template<>
int TPZCompElHDiv<TPZShapeCube>::ClassId() const
{
	return TPZHDIVCUBEID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapeCube>, TPZHDIVCUBEID>;

template<>
int TPZCompElHDiv<TPZShapeTetra>::ClassId() const
{
	return TPZHDIVTETRAID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapeTetra>, TPZHDIVTETRAID>;

template<>
int TPZCompElHDiv<TPZShapePrism>::ClassId() const
{
	return TPZHDIVPRISMID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapePrism>, TPZHDIVPRISMID>;

template<>
int TPZCompElHDiv<TPZShapePiram>::ClassId() const
{
	return TPZHDIVPYRAMID;
}

template class
TPZRestoreClass< TPZCompElHDiv<TPZShapePiram>, TPZHDIVPYRAMID>;


template class TPZCompElHDiv<TPZShapeTriang>;
template class TPZCompElHDiv<TPZShapePoint>;
template class TPZCompElHDiv<TPZShapeLinear>;
template class TPZCompElHDiv<TPZShapeQuad>;
template class TPZCompElHDiv<TPZShapeTetra>;
template class TPZCompElHDiv<TPZShapePrism>;
template class TPZCompElHDiv<TPZShapePiram>;
template class TPZCompElHDiv<TPZShapeCube>;


TPZCompEl * CreateHDivPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv<TPZShapePoint>(mesh,gel,index);
}


TPZCompEl * CreateHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDivBound2< TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeCube >(mesh,gel,index);
}

TPZCompEl * CreateHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapePrism>(mesh,gel,index);
}

TPZCompEl * CreateHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapePiram >(mesh,gel,index);
}

TPZCompEl * CreateHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
	return new TPZCompElHDiv< TPZShapeTetra >(mesh,gel,index);
}

