//
// Created by Michel Lesoinne on 3/28/18.
//

#ifndef FEM_VECTORREFERENCE_H
#define FEM_VECTORREFERENCE_H

#include <cstddef>

namespace FetiLib {

/// \brief A vector reference. Equivalent to a pointer with a length.
/// \details Vector Reference is used to pass any vector-type data to FetiLib.
/// It allows for runtime check of vector length. To call a method that requires
/// a VectorReference<double> for example, use {data, length} as arguments.
/// Use Scalar = const T for immutable data. For example VectorReference<const double>
template <typename Scalar>
class VectorReference {
public:
	/// \brief Constructor.
	VectorReference(Scalar *data, size_t length) : _data(data), length(length) {}
	/// \brief Access data at a given index.
	Scalar &operator[](int i) { return _data[i]; }
	Scalar operator[](int i) const { return _data[i]; }
	/// \brief Obtain a pointer to the data.
	Scalar *data() const { return _data; }
	/// \brief Size of the vector data.
	size_t size() const { return length; }

	const Scalar *begin() const { return _data; }
	const Scalar *end() const { return _data+length; }
private:
	Scalar *_data;
	size_t length;
};

} // namespace FetiLib

#endif //FEM_VECTORREFERENCE_H
