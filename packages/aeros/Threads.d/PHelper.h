#ifndef _PHELPER_H_
#define _PHELPER_H_

#include <memory>
#include <vector>
#include <Timers.d/GetTime.h>
#include <Threads.d/Paral.h>
#include <Timers.d/DistTimer.h>


template <typename FType>
class FunctorExecuter : public TaskDescr {
public:
	explicit FunctorExecuter(const FType &ft) : fctor(ft) {}

	FunctorExecuter(const FunctorExecuter &) = default;

	void runFor(int) override;
private:
	const FType &fctor; //!< The functor to execute.
};

template <typename FType>
void
FunctorExecuter<FType>::runFor(int i) {
	fctor(i);
}

template <typename FType>
auto makeExecuter(const FType &ftor) {
	return FunctorExecuter<FType>{ftor};
}


template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void execParal(int n, TA *target, void (TB::*f)(int, FArgs ...) const, Args &&...args)
{
	auto call =[&](int i) { (static_cast<const TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void execParal(int n, TA *target, void (TB::*f)(int, FArgs ...), Args &&...args)
{
	auto call =[&](int i) { (static_cast<TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void timedParal(DistTimer &timer, int n, TA *target, void (TB::*f)(int, FArgs ...), Args &&...args)
{
	double initTime = getTime();
	long initMem  = threadManager->memoryUsed();
	auto call =[&](int i) { (static_cast<TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execTimedParal(timer, n, &fe);
	timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
};


template <typename TA, typename TB, typename ... FArgs, typename ... Args, typename X = typename std::enable_if<sizeof...(FArgs) == sizeof...(Args)>::type>
void timedParal(DistTimer &timer, int n, TA *target, void (TB::*f)(int, FArgs ...) const, Args &&...args)
{
	auto call =[&](int i) { (static_cast<const TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	double initTime = getTime();
	long initMem  = threadManager->memoryUsed();
	threadManager->execTimedParal(timer, n, &fe);
	timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
};

namespace thread_details {

template<class A, class B, bool isIndexed = !std::is_convertible<B,A>::value>
struct DirectOrIndexed {
	inline static B subEval(B t, int) { return t; }
};

template<class A, class B>
struct DirectOrIndexed<A,B,true> {
	inline static A subEval(B t, int i) { return t[i]; }
};

}

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApply(int n, A * const *target, void(B::*fct)(Args ...),PassedArgs &&...pargs) {
	auto call =[&](int i) { (
		static_cast<B *>(target[i])->*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(std::forward<PassedArgs>(pargs), i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApply(int n, A * const *target, void(B::*fct)(Args ...) const,PassedArgs &&...pargs) {
	auto call =[&](int i) { (
		static_cast<B *>(target[i])->*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(std::forward<PassedArgs>(pargs), i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};


template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApply(const std::vector<A*> &target, void(B::*fct)(Args ...),PassedArgs &&...pargs) {
	auto call =[&](int i) { (
			static_cast<B *>(target[i])->*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(std::forward<PassedArgs>(pargs), i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(target.size(), &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApply(const std::vector<A*> &target, void(B::*fct)(Args ...) const,PassedArgs &&...pargs) {
	auto call =[&](int i) { (
			static_cast<B *>(target[i])->*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(std::forward<PassedArgs>(pargs), i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(target.size(), &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApply(const std::vector<std::unique_ptr<A>> &target, void(B::*fct)(Args ...),PassedArgs &&...pargs) {
	auto call =[&](int i) { (
		static_cast<B *>(target[i].get())->*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(std::forward<PassedArgs>(pargs), i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(target.size(), &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApply(const std::vector<std::unique_ptr<A>> &target, void(B::*fct)(Args ...) const,PassedArgs &&...pargs) {
	auto call =[&](int i) { (
		static_cast<const B *>(target[i].get())->*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(std::forward<PassedArgs>(pargs), i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(target.size(), &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, A *target, void(B::*fct)(Args ...),PassedArgs ...pargs) {
	auto call =[&](int i) {
		(static_cast<B &>(target[i]).*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(pargs, i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};


template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, A *target, void(B::*fct)(Args ...) const,PassedArgs ...pargs) {
	auto call =[&](int i) {
		(static_cast<B &>(target[i]).*fct)(thread_details::DirectOrIndexed<Args, PassedArgs>::subEval(pargs, i)...);
	};
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, A **target, void(B::*fct)(Args ...),PassedArgs ...pargs) {
	paralApply(n, target, fct, pargs...);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, A **target, void(B::*fct)(Args ...) const,PassedArgs ...pargs) {
	paralApply(n, target, fct, pargs...);
};


template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, std::vector<A*> &target, void(B::*fct)(Args ...),PassedArgs ...pargs) {
	paralApply(n, target.data(), fct, pargs...);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(int n, std::vector<A*> &target, void(B::*fct)(Args ...) const,PassedArgs ...pargs) {
	paralApply(n, target.data(), fct, pargs...);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(std::vector<A*> &target, void(B::*fct)(Args ...),PassedArgs ...pargs) {
	paralApply(static_cast<int>(target.size()), target.data(), fct, pargs...);
};

template <typename A, typename B, typename ... Args, typename ... PassedArgs>
void paralApplyToAll(std::vector<A*> &target, void(B::*fct)(Args ...) const,PassedArgs ...pargs) {
	paralApply(static_cast<int>(target.size()), target.data(), fct, pargs...);
};
#endif
