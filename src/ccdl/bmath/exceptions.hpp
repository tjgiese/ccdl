#ifndef _bmath_exceptions_hpp_
#define _bmath_exceptions_hpp_

#include <exception>
#include <string>

namespace ccdl
{
  class SingularMatrixException : public std::exception
  {
  public:
    SingularMatrixException( char const * msg ) : _what(msg) {};
    virtual ~SingularMatrixException() throw() {};
    virtual void rethrow() const { throw *this; };
    char const * what() const throw() { return _what.c_str(); };
  protected:
    std::string _what;
  };
}

#endif
