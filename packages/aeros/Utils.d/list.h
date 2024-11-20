#ifndef _LIST_H_
#define _LIST_H_

template<class T>
class Item {
public:
  Item( T Value);
  void setNext( Item *Next);
  T getValue() const;
  Item *getNext() const;
private:
  T value;
  Item *next;
};

template<class T>
class List {
public:
  List();
  ~List();
  void reset();
  void  insertItemAtEnd( T Value);
  void  insertItemAtStart( T Value);
  void  deleteItemAtStart();
  void  deleteItemAtEnd();
  Item<T> *getFirst() const;
  Item<T> *getLast() const;
  int getNumberOfItems() const;
private:
  int      numberOfItems;
  Item<T> *first;
  Item<T> *last;
};

template<class T>
class Iterator {
public:
  Iterator( const List<T> &List);
  void reset( const List<T> &List);
  Item<T> *getCurrent() const;
  bool more() const;
  void advanceCurrent();
private:
  Item<T> *current;
};

#ifdef _TEMPLATE_FIX_
#include <Utils.d/list.C>
#endif
#endif

