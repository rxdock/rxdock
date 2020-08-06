#ifndef SINGLETON_H
#define SINGLETON_H

namespace rxdock {

template <class singleton_type> class Singleton {
protected:
  enum action { DESTROY, CREATE, GET };
  static singleton_type *object(const action &a) {
    static singleton_type *obj = nullptr;
    switch (a) {
    case DESTROY:
      delete obj;
      break;
    case GET:
      return obj;
    case CREATE:
      return (obj = new singleton_type);
    }
    return nullptr;
  }

public:
  inline static singleton_type &instance() {
    if (object(GET) != nullptr) {
      return *object(GET);
    } else {
      return *object(CREATE);
    }
  }

  virtual ~Singleton() { object(DESTROY); }

protected:
  Singleton() {}
};

} // namespace rxdock

#endif
