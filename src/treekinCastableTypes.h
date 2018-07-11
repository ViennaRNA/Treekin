#ifndef _TREEKINCASTABLETYPES_
#define _TREEKINCASTABLETYPES_


#ifdef WITH_MPACK_GMP
# include <gmpxx.h>
#endif

#ifdef WITH_MPACK_QD
# include <qd/qd_real.h>
# include <qd/dd_real.h>
#endif

#ifdef WITH_MPACK_MPFR
# include <mpack/mpreal.h>
# include <mpack/mutils_mpfr.h>
#endif

#ifdef WITH_MPACK___FLOAT128
# include <mpack/mutils___float128.h>
#endif

namespace treekinCastableTypes {
#ifdef WITH_MPACK_QD

class qd_real_castable : public qd_real {
public:
  qd_real_castable() : qd_real()
  {}


  qd_real_castable(const qd_real & a) : qd_real(a)
  {}


  qd_real_castable(double a)
  {
    this->x[0] = a;
  }


  operator double() {
    return this->x[0];
  }
};

#endif

#ifdef WITH_MPACK_DD

class dd_real_castable : public dd_real {
public:
  dd_real_castable() : dd_real()
  {}


  dd_real_castable(const dd_real & a) : dd_real(a)
  {}


  dd_real_castable(double a)
  {
    this->x[0] = a;
  }


  operator double() {
    return this->x[0];
  }
};

#endif

class mpf_real_castable : public mpf_class {
public:

  mpf_real_castable() : mpf_class()
  {}


  mpf_real_castable(const mpf_class & a) : mpf_class(a)
  {}


  mpf_real_castable&
  operator=(const mpf_real_castable& v)
  {
    if (this != &v) {
      if (&v != NULL && v.__get_mp()->_mp_d != NULL) {
        if (this->__get_mp()->_mp_d != NULL)
          mpf_set(this->__get_mp(), v.__get_mp());
        else
          mpf_init_set(this->__get_mp(), v.__get_mp());
      }
    }

    return *this;
  }


  mpf_real_castable&
  operator=(const double& v)
  {
    if (this->__get_mp()->_mp_d != NULL)
      mpf_set_d(this->__get_mp(), v);
    else
      mpf_init_set_d(this->__get_mp(), v);

    return *this;
  }


  mpf_real_castable(double d)
  {
    if (this->__get_mp() == NULL)
      mpf_init2(this->__get_mp(), mpf_get_default_prec());

    mpf_set_d(this->__get_mp(), d);
  }


  operator double() {
    return this->get_d();
  }
};

#ifdef WITH_MPACK_MPFR

class mpreal_castable : public mpreal {
public:
  mpreal_castable() : mpreal()
  {}


  mpreal_castable(const double& a) : mpreal(a)
  {}


  mpreal_castable(const mpreal & a) : mpreal(a)
  {}


  operator int() const
  {
    return (int)mpfr_get_ui((mpfr_ptr)this, default_rnd);
  }
};

#endif

}

#ifdef WITH_MPACK_MPFR

namespace std {
static treekinCastableTypes::mpreal_castable
abs(treekinCastableTypes::mpreal_castable a)
{
  return abs((mpfr::mpreal)a);
}


static treekinCastableTypes::mpreal_castable
pow(treekinCastableTypes::mpreal_castable a,
    treekinCastableTypes::mpreal_castable b)
{
  return pow((mpfr::mpreal)a, (mpfr::mpreal)b);
}


static treekinCastableTypes::mpreal_castable
sqrt(treekinCastableTypes::mpreal_castable a)
{
  return sqrt((mpfr::mpreal)a);
}
}
#endif

#endif
