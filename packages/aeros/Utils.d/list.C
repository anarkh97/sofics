#ifndef _LIST_C_
#define _LIST_C_

#include "list.h"

template<class T>
Item<T>::Item( T Value)
{
  value = Value;
  next = 0;
}

template<class T>
void
Item<T>::setNext( Item *Next)
{
  next = Next;
}

template<class T>
T
Item<T>::getValue() const
{
  return value;
}

template<class T>
Item<T> *
Item<T>::getNext() const
{
  return next;
}

template<class T>
List<T>::List()
{
  numberOfItems = 0;
  first = 0;
  last = 0;
}

template<class T>
List<T>::~List()
{
  while( first != 0) {
    Item<T> *tmp = first;
	first = first->getNext();
	delete tmp;
	numberOfItems--;
  }
  last = 0;
}

template<class T>
void
List<T>::reset()
{
  while( first != 0) {
    Item<T> *tmp = first;
	first = first->getNext();
	delete tmp;
	numberOfItems--;
  }
  last = 0;
}

template<class T>
void
List<T>::insertItemAtEnd( T Value)
{
  Item<T> *item = new Item<T>( Value);
  if ( first == 0) {
    first = item;
	last = item;
  }
  else {
    last->setNext( item);
    last = item;
  }
  numberOfItems++;
}

template<class T>
void  
List<T>::insertItemAtStart( T Value)
{
  Item<T> *item = new Item<T>( Value);
  if ( first == 0) {
    first = item;
	last = item;
  }
  else {
    item->setNext( first);
    first = item;
  }
  numberOfItems++;
}

template<class T>
void  
List<T>::deleteItemAtStart()
{
  if ( first != 0) {
    Item<T> *firstCopy = first;
    first = first->getNext();
    delete firstCopy;
    numberOfItems--;
  }
  if ( numberOfItems == 0) {
    first = 0;
    last = 0;
  }
}

template<class T>
void  
List<T>::deleteItemAtEnd()
{
  if ( first != 0) {
    Item<T> *prev = first;
    while( ( prev->getNext() != last) && ( prev->getNext() != 0))
      prev = prev->getNext();
    if ( last == first) {
      delete last;
    }
    else {
      delete last;
      last = prev;
      last->setNext( 0);
    }
    numberOfItems--;
  }
  if ( numberOfItems == 0) {
    first = 0;
    last = 0;
  }
}

template<class T>
Item<T> *
List<T>::getFirst() const
{
  return first;
}

template<class T>
Item<T> *
List<T>::getLast() const
{
  return last;
}

template<class T>
int
List<T>::getNumberOfItems() const
{
  return numberOfItems;
}

template<class T>
Iterator<T>::Iterator( const List<T> &List)
{
  current = List.getFirst();
}

template<class T>
void
Iterator<T>::reset( const List<T> &List)
{
  current = List.getFirst();
}

template<class T>
Item<T> *
Iterator<T>::getCurrent() const
{
  return current;
}

template<class T>
bool
Iterator<T>::more() const
{
  if (current != 0) 
    return true;
  else
    return false;
}

template<class T>
void
Iterator<T>::advanceCurrent()
{
  current = current->getNext();
}

#endif
