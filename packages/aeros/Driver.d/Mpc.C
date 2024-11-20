#include <Driver.d/Mpc.h>

// real constructor
LMPCTerm::LMPCTerm(int _nnum, int _dofnum, double _coef)
{ 
  nnum = _nnum; 
  dofnum = _dofnum; 
  coef.r_value = _coef; 
  isComplex = false; 
}

// complex constructor
LMPCTerm::LMPCTerm(int _nnum, int _dofnum, double _rcoef, double _icoef)
{ 
  nnum = _nnum; 
  dofnum = _dofnum; 
  coef.c_value = DComplex(_rcoef, _icoef); 
  isComplex = true; 
}

// copy constructors
LMPCTerm::LMPCTerm(const LMPCTerm &t)
{
  isComplex = t.isComplex;
  nnum = t.nnum;
  dofnum = t.dofnum;
  if(t.isComplex) coef.c_value = t.coef.c_value;
  else coef.r_value = t.coef.r_value;
}

LMPCTerm::LMPCTerm(LMPCTerm &t, double weight)
{ 
  isComplex = t.isComplex;
  nnum = t.nnum;
  dofnum = t.dofnum;
  if(t.isComplex) coef.c_value = t.coef.c_value * weight;
  else coef.r_value = t.coef.r_value * weight;
}

LMPCTerm::LMPCTerm(const LMPCTerm &t, bool _isComplex)
{ 
  isComplex = _isComplex;
  nnum = t.nnum;
  dofnum = t.dofnum;
  if(isComplex) {
    if(t.isComplex) coef.c_value = t.coef.c_value;
    else coef.c_value = DComplex(t.coef.r_value, 0.0);
  }
  else {
    if(t.isComplex) coef.r_value = t.coef.c_value.real();
    else coef.r_value = t.coef.r_value;
  }
}

// default constructor
LMPCTerm::LMPCTerm(bool _isComplex)
{
  isComplex = _isComplex;
  nnum = dofnum = 0;
  if(isComplex) coef.c_value = DComplex(0.0, 0.0);
  else coef.r_value = 0.0;
}

// convert a real term to a complex term if it is not already complex
void 
LMPCTerm::makeComplex()
{
  if(!isComplex) {
    coef.c_value = DComplex(coef.r_value, 0.0);
    isComplex = true;
  }
}


// real constructor
LMPCons::LMPCons(int _lmpcnum, double _rhs, LMPCTerm *term0)
 : type(0), lagrangeMult(-1), penalty(0)
{
  isComplex = false;
  lmpcnum = _lmpcnum;
  original_rhs.r_value = rhs.r_value = _rhs;
  nterms = 0;
  m_type = mpc::Equality;
  m_source = mpc::Undefined;
  if(term0) addterm(term0);
}

// complex constructor
LMPCons::LMPCons(int _lmpcnum, double rrhs, double irhs, LMPCTerm *term0)
 : type(0), lagrangeMult(-1), penalty(0)
{
  isComplex = true;
  lmpcnum = _lmpcnum;
  original_rhs.c_value = rhs.c_value = DComplex(rrhs, irhs);
  nterms = 0;
  m_type = mpc::Equality;
  m_source = mpc::Undefined;
  if(term0) addterm(term0);
}

// convert rhs and all terms to complex if not already
void
LMPCons::makeComplex()
{
  if(!isComplex) {
    rhs.c_value = DComplex(rhs.r_value, 0.0);
    for(int i=0; i<nterms; ++i) terms[i].makeComplex();
    isComplex = true;
  }
}

// add term
void 
LMPCons::addterm(LMPCTerm *term)
{
  int i=0;
  // --- Verify if term already exists
  while((i < nterms) && ((terms[i].nnum != term->nnum) | (terms[i].dofnum != term->dofnum))) i++;
  // if term not already implied, add term
  if(i == nterms) {
    if(!isComplex && term->isComplex) makeComplex(); // new term is complex & lmpc was real
    if(isComplex && !term->isComplex) term->makeComplex();  // new term was real & lmpc is complex
    //terms[nterms++] = *term;
    terms.push_back(*term); 
    nterms++;
  }
  // if term already implied, add coefficients
  else {
     if(isComplex) {
       if(term->isComplex) terms[i].coef.c_value += term->coef.c_value;
       else terms[i].coef.c_value += DComplex(term->coef.r_value, 0.0);
     }
     else {
       if(term->isComplex) {  // need to convert all terms & rhs to complex
         makeComplex();
         terms[i].coef.c_value += term->coef.c_value;
       }
       else terms[i].coef.r_value += term->coef.r_value;
    }
  }
}

void 
LMPCons::print()
{
  if(!isComplex) {
    std::cerr << "lmpcnum = " << lmpcnum << ", rhs = " << rhs.r_value << ", nterms = " << nterms << std::endl;
    for(int i=0; i<nterms; ++i)
      std::cerr << "  term " << i+1 << ": node " << terms[i].nnum+1 << "  dof "
                << terms[i].dofnum << "  coef " << terms[i].coef.r_value << std::endl;
  }
  else {
    std::cerr << "lmpcnum = " << lmpcnum << ", rhs = " << rhs.c_value << ", nterms = " << nterms << std::endl;
    for(int i=0; i<nterms; ++i)
      std::cerr << "  term " << i+1 << ": node " << terms[i].nnum+1 << "  dof "
                << terms[i].dofnum << "  coef " << terms[i].coef.c_value << std::endl;
  }
}

// data access functions
template<>
double
LMPCons::getRhs<double>()
{
  if(isComplex) return rhs.c_value.real();
  else return rhs.r_value;
}

template<>
DComplex
LMPCons::getRhs<DComplex>()
{
  if(isComplex) return rhs.c_value;
  else return DComplex(rhs.r_value, 0.0);
}

template<>
GenLMPCTerm<double> 
LMPCons::getTerm<double>(int i)
{
  GenLMPCTerm<double> ret;
  ret.nnum = terms[i].nnum;
  ret.dofnum = terms[i].dofnum;
  if(isComplex) ret.coef = terms[i].coef.c_value.real();
  else ret.coef = terms[i].coef.r_value;
  return ret;
}

template<>
GenLMPCTerm<DComplex>
LMPCons::getTerm<DComplex>(int i)
{
  GenLMPCTerm<DComplex> ret;
  ret.nnum = terms[i].nnum;
  ret.dofnum = terms[i].dofnum;
  if(isComplex) ret.coef = terms[i].coef.c_value;
  else ret.coef = DComplex(terms[i].coef.r_value, 0.0);
  return ret;
}

void
LMPCons::normalize()
{
  if(isComplex) {
    std::cerr << " *** WARNING: LMPCons::normalize() not implemented for complex coefficients \n";
    return;
  }
  double cnorm = 0.0;
  for(int j=0; j<nterms; j++)
    cnorm += terms[j].coef.r_value * terms[j].coef.r_value;
  cnorm = sqrt(cnorm);
  for(int j=0; j<nterms; j++)
    terms[j].coef.r_value /= cnorm;
  rhs.r_value /= cnorm;
}

