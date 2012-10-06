/**
 * @file
 * @brief Contains the implementation of the TPZIntRuleT3D methods. 
 */

#include "tpzintrulet3d.h"
#include "pzerror.h"
#include "pzvec.h"

TPZIntRuleT3D::TPZIntRuleT3D(int order) {
	if(order < 0 && order > NRULESTETRAHEDRA_ORDER){
		PZError << "TPZGaussRule creation precision = " << order << " not available\n";
		order = NRULESTETRAHEDRA_ORDER;
	}
	ComputingSymmetricCubatureRule(order);
}

REAL TPZIntRuleT3D::W(int i) const {	
	if(i>=0 && i<fNumInt)
		return fWeight[i];
	else {
		PZError << "ERROR(TPZIntRuleT3D::w) Out of bounds!!\n";
		return 0.0;
	}
}

TPZIntRuleT3D::~TPZIntRuleT3D(){
	fLocationKsi.Resize(0);
	fLocationEta.Resize(0);
	fLocationZeta.Resize(0);
	fWeight.Resize(0);
	fNumInt = 0;
}

void TPZIntRuleT3D::Loc(int i, TPZVec<REAL> &Points) const {
	if(i>=0 && i<fNumInt){
		Points[0] = fLocationKsi[i];
		Points[1] = fLocationEta[i];
		Points[2] = fLocationZeta[i];
		return;
	}
	else {
		PZError << "ERROR(TPZIntRuleT::loc) Out of bounds!!\n";
	}
}

/** Symmetric Quadrature rule for tetrahedra from:
 Linbo Zhang, Tao Cui and Hui Liu, "A SET OF SYMMETRIC QUADRATURE RULES ON 
 TRIANGLES AND TETRAHEDRA", J. of Comput. Mathematics, Vol.27, No.1, 2009, 89–96.
 */

/**
 In the article is presented the following table, order versus number of points:
 Order:		1	2	3	4	5	6	7
 NPoints:	1	4	8	14	14	24	36
 
 Order:		8	9	10	11	12	13	14
 NPoints:	46	61	81	109	140	171	236
 */

#ifdef Length
#undef Length
#endif
#define Length(wts)	(sizeof(wts) / (sizeof(wts[0])))

#define Perm4(a)	a,a,a,a
#define Dup4(w)		((1.0L/6.0L)*w)

#define Perm31(a)	a,a,a,(1.L-(3.L*(a))), \
a,a,(1.L-(3.L*(a))),a, \
a,(1.L-(3.L*(a))),a,a, \
(1.L-(3.L*(a))),a,a,a
#define Dup31(w)	Dup4(w),Dup4(w),Dup4(w),Dup4(w)

#define Perm22(a)	a,a,(0.5L-(a)),(0.5L-(a)), \
a,(0.5L-(a)),a,(0.5L-(a)), \
a,(0.5L-(a)),(0.5L-(a)),a, \
(0.5L-(a)),a,(0.5L-(a)),a, \
(0.5L-(a)),a,a,(0.5L-(a)), \
(0.5L-(a)),(0.5L-(a)),a,a
#define Dup22(w)	Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w)

#define Perm211(a,b)	a,a,b,(1.0L-(a)-(a)-(b)), \
a,a,(1.0L-(a)-(a)-(b)),b, \
a,b,a,(1.0L-(a)-(a)-(b)), \
a,b,(1.0L-(a)-(a)-(b)),a, \
a,(1.0L-(a)-(a)-(b)),a,b, \
a,(1.0L-(a)-(a)-(b)),b,a, \
b,a,a,(1.0L-(a)-(a)-(b)), \
b,a,(1.0L-(a)-(a)-(b)),a, \
b,(1.0L-(a)-(a)-(b)),a,a, \
(1.0L-(a)-(a)-(b)),a,a,b, \
(1.0L-(a)-(a)-(b)),a,b,a, \
(1.0L-(a)-(a)-(b)),b,a,a
#define Dup211(w)	Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w), \
	Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w)

#define Perm0111(p,a,b,c) p,a,b,c, p,a,c,b, p,b,a,c, p,b,c,a, p,c,a,b, p,c,b,a
#define Perm1111(a,b,c) Perm0111(a,b,c,(1.0L-(a)-(b)-(c))), \
Perm0111(b,a,c,(1.0L-(a)-(b)-(c))), \
Perm0111(c,a,b,(1.0L-(a)-(b)-(c))), \
Perm0111((1.0L-(a)-(b)-(c)),a,b,c)

// #define Dup111(w)	((1.L/6.L)*w), ((1.L/6.L)*w), ((1.L/6.L)*w), ((1.L/6.L)*w), ((1.L/6.L)*w), ((1.L/6.L)*w)

#define Dup1111(w)	Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w), \
	Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w), \
	Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w),Dup4(w)

long double QUAD_TETRAH_P1_wts[] = {
    Dup4(1.L)
};
long double QUAD_TETRAH_P1_pts[Length(QUAD_TETRAH_P1_wts) * 4] = {
    Perm4(.25L)
};

long double QUAD_TETRAH_P2_wts[] = {
    Dup31(.25L)
};
long double QUAD_TETRAH_P2_pts[Length(QUAD_TETRAH_P2_wts) * 4] = {
    /* (5. - sqrt(5.)) / 20, (5. + 3 * sqrt(5.)) / 20 */
    Perm31(.13819660112501051517954131656343619L)
};

long double QUAD_TETRAH_P3_wts[] = {
    /* 1 / 8 + sqrt((1715161837 - 406006699 * sqrt(17)) / 23101) / 3120 */
    Dup31(.13852796651186214232361769837564129L),
    /* 1 / 8 - sqrt((1715161837 - 406006699 * sqrt(17)) / 23101) / 3120 */
    Dup31(.11147203348813785767638230162435871L)
};
long double QUAD_TETRAH_P3_pts[Length(QUAD_TETRAH_P3_wts) * 4] = {
    /* (55 - 3 * sqrt(17) + sqrt(1022 - 134 * sqrt(17))) / 196 */
    Perm31(.32805469671142664733580581998119743L),
    /* (55 - 3 * sqrt(17) - sqrt(1022 - 134 * sqrt(17))) / 196 */
    Perm31(.10695227393293068277170204157061650L)
};

long double QUAD_TETRAH_P4_wts[] = {
    Dup31(.07349304311636194934358694586367885L),
    Dup31(.11268792571801585036501492847638892L),
    Dup22(.04254602077708146686093208377328816L)
};
long double QUAD_TETRAH_P4_pts[Length(QUAD_TETRAH_P4_wts) * 4] = {
    Perm31(.09273525031089122628655892066032137L),
    Perm31(.31088591926330060975814749494040332L),
    Perm22(.04550370412564965000000000000000000L)
};

long double QUAD_TETRAH_P5_wts[] = {
    Dup31(.11268792571801585079918565233328633L),
    Dup31(.07349304311636194954371020548632750L),
    Dup22(.04254602077708146643806942812025744L)
	
};
long double QUAD_TETRAH_P5_pts[Length(QUAD_TETRAH_P5_wts) * 4] = {
    Perm31(.31088591926330060979734573376345783L),
    Perm31(.09273525031089122640232391373703061L),
    Perm22(.04550370412564964949188052627933943L)
};

long double QUAD_TETRAH_P6_wts[] = {
    Dup31(.03992275025816749209969062755747998L),
    Dup31(.01007721105532064294801323744593686L),
    Dup31(.05535718154365472209515327785372602L),
    Dup211(27./560.L)
};
long double QUAD_TETRAH_P6_pts[Length(QUAD_TETRAH_P6_wts) * 4] = {
    Perm31(.21460287125915202928883921938628499L),
    Perm31(.04067395853461135311557944895641006L),
    Perm31(.32233789014227551034399447076249213L),
    /* (3 - sqrt(5)) / 12, (5 + sqrt(5)) / 12, (1 + sqrt(5)) / 12 */
    Perm211(.06366100187501752529923552760572698L,
			.60300566479164914136743113906093969L)
};

long double QUAD_TETRAH_P7_wts[] = {
    Dup4(.09548528946413084886057843611722638L),
    Dup31(.04232958120996702907628617079854674L),
    Dup22(.03189692783285757993427482408294246L),
    Dup211(.03720713072833462136961556119148112L),
    Dup211(.00811077082990334156610343349109654L)
};
long double QUAD_TETRAH_P7_pts[Length(QUAD_TETRAH_P7_wts) * 4] = {
    Perm4(.25),
    Perm31(.31570114977820279942342999959331149L),
    Perm22(.05048982259839636876305382298656247L),
    Perm211(.18883383102600104773643110385458576L,
			.57517163758700002348324157702230752L),
    Perm211(.02126547254148324598883610149981994L,
			.81083024109854856111810537984823239L)
};

long double QUAD_TETRAH_P8_wts[] = {
    Dup31(.00639714777990232132145142033517302L),
    Dup31(.04019044802096617248816115847981783L),
    Dup31(.02430797550477032117486910877192260L),
    Dup31(.05485889241369744046692412399039144L),
    Dup22(.03571961223409918246495096899661762L),
    Dup211(.00718319069785253940945110521980376L),
    Dup211(.01637218194531911754093813975611913L)
};
long double QUAD_TETRAH_P8_pts[Length(QUAD_TETRAH_P8_wts) * 4] = {
    Perm31(.03967542307038990126507132953938949L),
    Perm31(.31448780069809631378416056269714830L),
    Perm31(.10198669306270330000000000000000000L),
    Perm31(.18420369694919151227594641734890918L),
    Perm22(.06343628775453989240514123870189827L),
    Perm211(.02169016206772800480266248262493018L,
			.71993192203946593588943495335273478L),
    Perm211(.20448008063679571424133557487274534L,
			.58057719012880922417539817139062041L)
};

long double QUAD_TETRAH_P9_wts[] = {
    Dup4(.05642669317950620658871504327612541L),
    Dup31(.00334109507471348040299974430471765L),
    Dup31(.03011375476877376390731423843157491L),
    Dup31(.00649096092006153463576211689456861L),
    Dup211(.00980928586825458643196874259255500L),
    Dup211(.02811915382336547255163261742529262L),
    Dup211(.00789458690833150076834149200960885L),
    Dup211(.01949281204723999671697219448924602L)
};
long double QUAD_TETRAH_P9_pts[Length(QUAD_TETRAH_P9_wts) * 4] = {
    Perm4(.25),
    Perm31(.03402217700104486646540370887876764L),
    Perm31(.32277033353380052539137668325496398L),
    Perm31(.06045707742577493000000000000000000L),
    Perm211(.45536299094720821180030815044164301L,
			.00568317736533017990610016014574474L),
    Perm211(.11950225539382580097797370469611438L,
			.46311683247848994097622449365772955L),
    Perm211(.02802195578340115815505750665412373L,
			.72520607683986748873856595428480993L),
    Perm211(.17483303201157461578532464597224522L,
			.61668257178125640457068309097954073L)
};

long double QUAD_TETRAH_P10_wts[] = {
    Dup4(.04739977355602073838473882117805110L),
    Dup31(.02693705999226869980276416100488208L),
    Dup31(.00986915971679338323455773543017308L),
    Dup211(.00036194434433925362423987838480851L),
    Dup211(.01013587167975579278851647011501678L),
    Dup211(.01139388122019523162362093488071434L),
    Dup211(.00657614727703590416745574020045070L),
    Dup211(.02573973198045607127903601225965471L),
    Dup211(.01290703579886199063929543024949899L)
};
long double QUAD_TETRAH_P10_pts[Length(QUAD_TETRAH_P10_wts) * 4] = {
    Perm4(.25),
    Perm31(.31225006869518864772980831868682746L),
    Perm31(.11430965385734615058737119765365045L),
    Perm211(.00613800882479074784759371324841535L,
			.94298876734520486619763058691825076L),
    Perm211(.03277946821644267077472102033232419L,
			.34018479408710763278898792494967132L),
    Perm211(.41043073921896549428789784425151169L,
			.16548602561961105160449012444452641L),
    Perm211(.03248528156482304783551493997842620L,
			.13385215221200951309782843596456662L),
    Perm211(.12105018114558942599389500159505053L,
			.47719037990428035054410640829690722L),
    Perm211(.17497934218393902428494922652831040L,
			.62807184547536601069327607221790967L)
};

long double QUAD_TETRAH_P11_wts[] = {
    Dup4(.03943210802865886350733033449120443L),
    Dup31(.01566212622727911315008856276876506L),
    Dup31(.00333217237490140814440923615401491L),
    Dup31(.01402607740748974743749136099769235L),
    Dup211(.00108590752933246630682209837723547L),
    Dup211(.02023596043066317891111657316540838L),
    Dup211(.01179021487212586353684938046770181L),
    Dup211(.00769031498252129590113157802073890L),
    Dup211(.00443730570345920390473072602143959L),
    Dup211(.01142954846718404041077055259859402L),
    Dup1111(.00618564017121781141281925508389534L)
};
long double QUAD_TETRAH_P11_pts[Length(QUAD_TETRAH_P11_wts) * 4] = {
    Perm4(.25L),
    Perm31(.12149136777653379449770230990807224L),
    Perm31(.03231625915107289635395445208958103L),
    Perm31(.32492614978860679781284190241442197L),
    Perm211(.00414835697166001200000000000001000L,
			.59826599679018635020545384277617780L),
    Perm211(.22462461067637714141447515116498644L,
			.47366228783234957140836966920205236L),
    Perm211(.05190508777256569674422721644265892L,
			.56314477790827989873710197630305713L),
    Perm211(.13493013121624020422375917234299303L,
			.70835883078581895385699500512712996L),
    Perm211(.02519119210825247292005118506530550L,
			.78371950734007737543057403429990901L),
    Perm211(.36531877978173361396933198009886720L,
			.13460390831686580000000000000001000L),
    Perm1111(.52290753950993847296521692758602923L,
			 .14075363054369590184253913949127849L,
			 .00976243819645261550829228038997777L)
};

long double QUAD_TETRAH_P12_wts[] = {
    Dup31(.01276763770097074150203778596512505L),
    Dup31(.01612110423790926821858154489575762L),
    Dup31(.00037161269857844220004255818986081L),
    Dup31(.01971744178668545763955330903818868L),
    Dup31(.00257139093086271836218234759448548L),
    Dup22(.00381724787051057590575318412783326L),
    Dup22(.01208722707766311317860318419314605L),
    Dup211(.00310586115843473343431688149929620L),
    Dup211(.00545953133647103066912742126769441L),
    Dup211(.00214289974849699750666852093655947L),
    Dup211(.00552467146725782962244930098165075L),
    Dup211(.00853695669449918042985177836672201L),
    Dup1111(.01151017784832330697333644123403294L),
    Dup1111(.00520387865288561360396792421252454L)
};
long double QUAD_TETRAH_P12_pts[Length(QUAD_TETRAH_P12_wts) * 4] = {
    Perm31(.11529974435148014530455720738915911L),
    Perm31(.20233628224059090000000000000001000L),
    Perm31(.01171759795761995151247906754831398L),
    Perm31(.31330644136780106727760279964458934L),
    Perm31(.25000573011558370000000000000001000L),
    Perm22(.02099547435075800669020182527059018L),
    Perm22(.15177401824745010000000000000001000L),
    Perm211(.02441977874343536478314000904761661L,
			.84832928469787285064520886743481574L),
    Perm211(.25620709853201830896382010708562210L,
			.48248737387384884780289289672973542L),
    Perm211(.01679032097960299061471796028857942L,
			.69477194236575592695949850988417719L),
    Perm211(.12616082113987204239970703846895919L,
			.72541048930294811897485950521263380L),
    Perm211(.43143517452637984721670695066371957L,
			.11272193989285241520959977211007542L),
    Perm1111(.50167006246250569747515507168476130L,
			 .27247180286952239178351046753060445L,
			 .07207432880729891465015948456335820L),
    Perm1111(.26164485453781874566945505006396799L,
			 .08629229194706173191742351944352488L,
			 .02056541065587613830062489762710900L)
};

long double QUAD_TETRAH_P13_wts[] = {
    Dup4(.01501368777308314675062970631615983L),
    Dup31(.01822520928017342532379068941490097L),
    Dup31(.00700610921774146424038518693926311L),
    Dup22(.01642354974394954829540573107905531L),
    Dup22(.00512061009636059707262596949702171L),
    Dup22(.01119669865290491634382032086351956L),
    Dup211(.01561914973337995400953811302431969L),
    Dup211(.00248442301331647441904056776338473L),
    Dup211(.00163859853481823893844525309440751L),
    Dup211(.00590303044012492197171914655535865L),
    Dup211(.01102208245821805240445097989201525L),
    Dup211(.00040645183996417822585155512755848L),
    Dup1111(.00268796997296854209745781926651729L),
    Dup1111(.00197950480552671190531894675510740L),
    Dup1111(.00544631918142579120943187040108667L)
};
long double QUAD_TETRAH_P13_pts[Length(QUAD_TETRAH_P13_wts) * 4] = {
    Perm4(.25),
    Perm31(.15521609351908950314115784335704739L),
    Perm31(.33012266333967360024433192595196779L),
    Perm22(.16680640389386249928937782601144234L),
    Perm22(.02492378854777361779701400374860089L),
    Perm22(.09719762991575100143072243716240818L),
    Perm211(.24785929015736256692746910620827934L,
			.43365324235685144718726061434767377L),
    Perm211(.02223159608186700290879521860892929L,
			.83690032040373400514509486595698594L),
    Perm211(.10727869331305341049150459639584801L,
			.77498030597500180756587877274179289L),
    Perm211(.19817684388398981142331840582142759L,
			.58756930578220530259172017903595920L),
    Perm211(.06917924347737931647732534347465502L,
			.60420006666006644707935264871115302L),
    Perm211(.02311471947193316000000000000001000L,
			.93087579279244424864920228882888307L),
    Perm1111(.11788928751019608922290117470644250L,
			 .11651536422540720000000000000001000L,
			 .04202400112551542095676634303719997L),
    Perm1111(.67703279860228426355032221326746594L,
			 .04616537602461971083458041122176081L,
			 .00084434031890503975729899692135905L),
    Perm1111(.48489008867363312201080094154790828L,
			 .35888294295520201572423646909421086L,
			 .13818283491762872996955080907912355L)
};

long double QUAD_TETRAH_P14_wts[] = {
    Dup31(.00406511366527076704362088368356360L),
    Dup31(.00221453853344557814375995695000715L),
    Dup31(.00581343826788845054953733388214554L),
    Dup31(.01962554338583572159756233339617148L),
    Dup31(.00038757379059082143645387212483937L),
    Dup211(.01164297197217703698552134010055516L),
    Dup211(.00528904298828171313177368830528561L),
    Dup211(.00183108541636005593766978234880692L),
    Dup211(.00824964737721464520674496691736603L),
    Dup1111(.00300992453470824513768887482089866L),
    Dup1111(.00080471656173675346362618087603116L),
    Dup1111(.00298504125884930711876556928839215L),
    Dup1111(.00568960024187607669633614778119730L),
    Dup1111(.00415908658785457156700139801826135L),
    Dup1111(.00072823892045727243561364297456536L),
    Dup1111(.00543265007699582482162423406519264L)
};
long double QUAD_TETRAH_P14_pts[Length(QUAD_TETRAH_P14_wts) * 4] = {
    Perm31(.32725336252384856390930966926852893L),
    Perm31(.04476130446668508088379420964788419L),
    Perm31(.08614033110243635365372087402988575L),
    Perm31(.20876264250043229682653570839761758L),
    Perm31(.01410497380292096006358791521029282L),
    Perm211(.10216532418077681234766925269825839L,
			.57394636759433382028140028934601068L),
    Perm211(.40757005166001071572132956513017833L,
			.09222787013902013000000000000000000L),
    Perm211(.01566400074028035855575867095780840L,
			.70128109595894403271399676732084261L),
    Perm211(.22549635625250290537807241542011034L,
			.47690639744208871158605833541070112L),
    Perm1111(.39059842812814580000000000000000000L,
			 .20135905441239221681230773272350923L,
			 .01611228807103002985780269315483708L),
    Perm1111(.10613506799890214555561390298480794L,
			 .03273581868172692849440040779126601L,
			 .00359790765372716669079715233859245L),
    Perm1111(.56363837316977438968968166306485017L,
			 .23029207223006574545025268741356515L,
			 .19071993417435518627124877906378985L),
    Perm1111(.36762550953258608440922067759911669L,
			 .20788513802300449507171021252507348L,
			 .33121048851934490000000000000000000L),
    Perm1111(.71923236898172952950234018407969909L,
			 .17632791180193297621579930336369727L,
			 .02076023625713100907549734406116442L),
    Perm1111(.52782499521529872984092400758172763L,
			 .43728908922034181655262387608419181L,
			 .00922016518566419494631775549492202L),
    Perm1111(.54836745449481907289949105056077457L,
			 .34478155061716412287036718709203314L,
			 .08672172833222153946294387400858277L)
};

void TPZIntRuleT3D::ComputingSymmetricCubatureRule(int order) {
	if(order > 14) order = 14; 
	int NRGAUPO[15] = { 1, 1, 4, 8, 14, 14, 24, 36, 46, 61, 81, 109, 140, 171, 236};
	fNumInt = NRGAUPO[order];
	fLocationKsi.Resize(fNumInt,0.0L);
	fLocationEta.Resize(fNumInt,0.0L);
	fLocationZeta.Resize(fNumInt,0.0L);
	fWeight.Resize(fNumInt,0.0L);
	
	switch(order) {
		case 0:
		case 1:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P1_pts,QUAD_TETRAH_P1_wts);
			break;
		case 2:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P2_pts,QUAD_TETRAH_P2_wts);
			break;
		case 3:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P3_pts,QUAD_TETRAH_P3_wts);
			break;
		case 4:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P4_pts,QUAD_TETRAH_P4_wts);
			break;
		case 5:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P5_pts,QUAD_TETRAH_P5_wts);
			break;
		case 6:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P6_pts,QUAD_TETRAH_P6_wts);
			break;
		case 7:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P7_pts,QUAD_TETRAH_P7_wts);
			break;
		case 8:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P8_pts,QUAD_TETRAH_P8_wts);
			break;
		case 9:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P9_pts,QUAD_TETRAH_P9_wts);
			break;
		case 10:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P10_pts,QUAD_TETRAH_P10_wts);
			break;
		case 11:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P11_pts,QUAD_TETRAH_P11_wts);
			break;
		case 12:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P12_pts,QUAD_TETRAH_P12_wts);
			break;
		case 13:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P13_pts,QUAD_TETRAH_P13_wts);
			break;
		case 14:
			TransformBarycentricCoordInCartesianCoord(QUAD_TETRAH_P14_pts,QUAD_TETRAH_P14_wts);
			break;
		default:
			PZError << "TPZIntRuleT3D not implemented by order " << order << std::endl;
	}
}

void TPZIntRuleT3D::TransformBarycentricCoordInCartesianCoord(long double baryvec[],long double weightvec[]) {
	for(int i=0;i<fNumInt;i++) {
		fWeight[i] = weightvec[i];
		fLocationKsi[i] = baryvec[4*i+1];
		fLocationEta[i] = baryvec[4*i+2];
		fLocationZeta[i] = baryvec[4*i+3];
	}
}


