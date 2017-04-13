
template<typename MyGrid,typename MyFFT,typename AO,typename BO> void ABBAstepT(double tau, MyFFT& DFTval,MyGrid&  uIn,MyGrid&  uOut,MyGrid& potential,AO Aop,BO Bop,double alpha  = -1.0)
{	
	//AB step
	MyGrid  AUX=uIn; 
	AUX.setToValue(0.0);
	(*Aop)(tau/2.0,DFTval,uIn,AUX,alpha);
    MyGrid AUX2=uIn; 
	AUX2.setToValue(0.0);
	(*Bop)(tau,AUX, AUX2,potential);
	//BA step
	Aop(tau/2.0,DFTval,AUX2,uOut,alpha);

};
//templated ABBA step for A step, B step. A requires a DFT
template<typename MyGrid,typename MyFFT,typename AO,typename BO> void ABBAstepwithProjT(double tau, MyFFT& DFTval,MyGrid&  uIn,MyGrid&  uOut,MyGrid& potential,AO Aop,BO Bop)
{		ABBAstepT(tau,DFTval,uIn,uOut,potential,Aop,Bop);
	uOut/=uOut.norm2();
};

template<typename MyGrid, typename MyFFT,typename AO,typename BO> void RichardsonT(double tau, MyFFT& DFTval,MyGrid&  uIn,MyGrid&  uOut,MyGrid& potential,AO Aop,BO Bop,double alpha  = -1.0)
{	MyGrid AUX=uIn;
	AUX.setToValue(0.0);
	MyGrid AUX2=uIn; 
	AUX.setToValue(0.0);
	uOut.setToValue(0.0);
	double alph=-1.0/7.0;
	double bet=8.0/7.0;
	ABBAstepT(tau,DFTval,uIn,AUX,potential,Aop,Bop);
	//two half steps
	ABBAstepT(tau/2.0,DFTval,uIn,AUX2,potential,Aop,Bop);
	ABBAstepT(tau/2.0,DFTval,AUX2,uOut,potential,Aop,Bop);
	uOut*=bet;
	AUX*=alph;
	uOut+=AUX;
	uOut/=uOut.norm2();	
}

template<typename MyGrid, typename MyFFT,typename AO,typename BO> void RichardsonRecNOprojT(double tau, MyFFT& DFTval,MyGrid&  uIn,MyGrid&  uOut,MyGrid& potential,AO Aop,BO Bop,int order,double alpha  = -1.0){
	MyGrid AUX=uIn;
	AUX.setToValue(0.0);
	MyGrid AUX2=uIn; 
	AUX.setToValue(0.0);
	uOut.setToValue(0.0);
	double alph=-1.0/(pow(2.0,order+1)-1);
	double bet=1+alpha;
	if(order<=1)
		ABBAstepT(tau,DFTval,uIn,uOut,potential,Aop,Bop);
	else{//full step with lower order
		RichardsonRecNOprojT(tau, DFTval, uIn,AUX,potential,Aop,Bop,order-1,alpha);
		//two half steps
		RichardsonRecNOprojT(tau/2.0, DFTval, uIn,AUX2,potential,Aop,Bop,order-1,alpha);
		RichardsonRecNOprojT(tau/2.0, DFTval,AUX2,uOut,potential,Aop,Bop,order-1,alpha);
		uOut*=bet;
	AUX*=alph;
	uOut+=AUX;
		}	
}
template<typename MyGrid, typename MyFFT,typename AO,typename BO> void RichardsonRecANDprojT(double tau, MyFFT& DFTval,MyGrid&  uIn,MyGrid&  uOut,MyGrid& potential,AO Aop,BO Bop,int order,double alpha  = -1.0){
			RichardsonRecNOprojT(tau, DFTval, uIn,uOut,potential,Aop,Bop,order,alpha);
	
			uOut/=uOut.norm2();		
	}