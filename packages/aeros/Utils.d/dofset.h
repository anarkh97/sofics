#ifndef _DOFSET_H_
#define _DOFSET_H_

#include <iostream>
#include <vector>

class DofSet {
protected:
	int flags;
public:
	static const int max_known_dof = 24;
	static const int max_known_nonL_dof = 9;
	static int
		Xdisp,
		Ydisp,
		Zdisp,
		XYZdisp,
		Xrot,
		Yrot,
		Zrot,
		XYZrot,
		Temp,
		Helm,
		Contact,
		IntPress,
		Potential,
		LagrangeE,
		LagrangeI;
	static int nonL_dof;

	static const int DispAndRot = 0x3f;
	static DofSet nullDofset;

	// constructors
	DofSet() { flags = 0; }

	DofSet(int t) { flags = t; }

	DofSet &operator|=(const DofSet &ds) {
		flags |= ds.flags;
		return *this;
	}

	bool operator==(const DofSet &ds) const { return flags == ds.flags; }

	// mark marks given dofs as being used
	void mark(int dofs) { flags |= dofs; }

	void unmark(int dofs) { //fprintf(stderr,"In unmark: ");
		//fprintf(stderr,"flags = %d; ", flags);
		//fprintf(stderr,"dofs = %d\n", dofs);
		flags = flags ^ (flags & dofs);
	}

	void unmarkAll() { flags = 0; }

	// tests if all the given dofs are present
	int test(int dofs) { return (flags & dofs) == dofs; }

	/** Check presence of at least one dof among a given set.
	 *
	 * @param dofs DOFs whose presence should be checked.
	 * @return Whether any of the dofs are present.
	 */
	bool contains(int dofs) const {
		return (flags & dofs) != 0;
	}

/*
    bool containsAllDisp(int dim) {
      if(dim >=3) return contains(XYZdisp);
      else if(dim == 2) return contains(Xdisp & Ydisp);
      else if(dim == 1) return contains(Xdisp);
      else return true;
    }
*/
	bool containsAllDisp(int dim) const {
		switch (dim) {
			case 1 :
				return contains(Xdisp);
				break;
			case 2:
				return (contains(Xdisp) && contains(Ydisp));
				break;
			case 3:
				return (contains(Xdisp) && contains(Ydisp) && contains(Zdisp));
				break;
			default:
				std::cerr << " *** WARNING: in DofSet::containsAllDisp(), dim = " << dim << std::endl;
				return false;
				break;
		}
	}

	bool containsAnyRot() const {
		return contains(Xrot | Yrot | Zrot);
	}

	/// \brief Count the number dofs marked in this node.
	int count() const;

	int number(DofSet, int *) const; // This routine complements the previous
	// one, numbering all the DOFs in the first argument

	// Locates a dof in the current set. Only works for single dofs
	// i.e. result is undefined for composites such as XYZrot and XYZdisp
	int locate(int dof) const;

	int list() const { return flags; }

	DofSet operator&(DofSet d) const { return {flags & d.flags}; }

	DofSet operator|(DofSet d) const { return {flags | d.flags}; }

	DofSet operator^(DofSet d) const { return {flags ^ d.flags}; }

	DofSet &operator&=(const DofSet &d) {
		flags &= d.flags;
		return *this;
	}
	void print(char *msg = 0) const;
};

class Elemset;

class Element;


class EqNumberer {
protected:
	int numnodes{};                 //!< \brief total number of nodes (Now redundant with vector sizes)
	std::vector<int> node_offset;   //!< \brief First DOF # associated with a node.
	std::vector<int> node_num_dofs; //!< \brief Number of dofs associated with a node.
	std::vector<int> renummap;      //!< \brief Renumbering mapping of the nodes.

public:
	EqNumberer() {}

	virtual ~EqNumberer() {}

	// Return the total number of degrees of freedom
	int size() const { return (node_offset.size() != 0) ? node_offset[numnodes] : 0; }

	// Return the number of nodes
	int numNodes() const { return numnodes; }

	// Return the first dof of a node
	int firstdof(int node) const { return (node_num_dofs[node] > 0) ? node_offset[node] : -1; }

	// Return the weight of a node (its number of dofs)
	int weight(int n) const { return node_num_dofs[n]; }

	int *allWeights() { return node_num_dofs.data(); }

	int *allOffsets() { return node_offset.data(); }

	const int *renumPtr() const { return renummap.data(); }
	// MLX Get rid of this one.
	int *renumPtr() { return renummap.data(); }

	void print();
};

/** \brief Set of DOFs that appear with the nodes of a problem, associating a unique number to each.
 * \details The DofSetArray keeps track of all the DOFs that are present in each node.
 * The construction does not follow a RAII principle.
 * Instead there is a construction phase and then a use phase. It would be good to get to a RAII approach.
 */
class DofSetArray : public EqNumberer {
protected:
	std::vector<DofSet> dofs;
	std::vector<int> rowcolnum;
	std::vector<int> invrowcol;
	std::vector<int> dofType; // 0 = translational, 1 = rotational

	DofSetArray() = default;

protected:
	void makeOffset();

	void makeModifiedOffset();

public:
	DofSetArray(int nnode, int *dofsPerNode, int *renumtable); // for DEC
	DofSetArray(int nnodes, Elemset &elearray, int *renumtable = 0, int myMap = 0);

	DofSetArray(int nnodes, int *renumtable = 0, int myMap = 0);

	DofSetArray(DofSetArray &&dofSetArray) = default;

	explicit DofSetArray(Element *ele);

	~DofSetArray() override;

	// locate a dof for a given node
	int locate(int node, int dof) const;

	int number(int node, DofSet, int *) const;

	/// Mark dofs for a node
	void mark(int node, int dof);

	/** \brief Mark DOFs used at a node. */
	void mark(int node, DofSet ds) { mark(node, ds.list()); }

	void mark(const int *node, int numNode, int dof);

	// Return the DofSet of a node
	DofSet &operator[](int i) { return dofs[i]; }

	const DofSet &operator[](int i) const { return dofs[i]; }

	/// \brief Get the unconstrained index for a full set dof index or -1 if it is constrained.
	int getRCN(int dof) const { return rowcolnum[dof]; }
	/// \brief Get the constrained index for a full set dof index or -1 if it is unconstrained.
	int invRCN(int dof) const { return invrowcol[dof]; }
	/// \brief Get direct access to the vector mapping from full to unconstrained
	auto &getUnconstrNum() const { return rowcolnum; }
	/// \brief Get direct access to the vector mapping from full to constrained
	auto &getConstrndNum() const { return invrowcol; }

	int *makeDofTypeArray(); // creates and returns dofType array

	void setWeight(int n, int w);

	int getWeight(int n) const;

	void finish() { makeModifiedOffset(); }

	friend class ConstrainedDSA;

	void clean_up();
};

template<typename VecType>
void
zeroRotDofs(const DofSetArray &dsa, VecType &vec) {
	static const int ROT_DOFS[] = {DofSet::Xrot, DofSet::Yrot, DofSet::Zrot};
	static const int ROT_DOFS_SIZE = sizeof(ROT_DOFS) / sizeof(ROT_DOFS[0]);

	const int nodeCount = dsa.numNodes();
	for (int iNode = 0; iNode < nodeCount; ++iNode) {
		for (const int *dofType = ROT_DOFS; dofType != ROT_DOFS + ROT_DOFS_SIZE; ++dofType) {
			const int dofLoc = dsa.locate(iNode, *dofType);
			if (dofLoc >= 0) {
				vec[dofLoc] = 0.0;
			}
		}
	}
}

class BCond;

class ComplexBCond;

/** \brief DOF set where some DOFs have been constrained --i.e. removed from the active set.
 *
 */
class ConstrainedDSA : public DofSetArray {
public:
	ConstrainedDSA(const DofSetArray &dsa,
	               int nBC, const BCond *boundaryCond,
	               int nComplexBC = 0, const ComplexBCond *complexBC = nullptr);

	ConstrainedDSA(const DofSetArray &dsa, const DofSetArray &c_dsa, int nbc, const BCond *bcd, int info);

	ConstrainedDSA(const DofSetArray &dsa, int nbc, const BCond *bcd, const int *bc);

	ConstrainedDSA(const DofSetArray &dsa, int nbc, const ComplexBCond *bcd,
	               const int *bc);

	ConstrainedDSA(const DofSetArray &dsa, const BCond *bcdata, int nbc,
	               const ComplexBCond *bcd, int nbcd, const int *bc);

	ConstrainedDSA(const DofSetArray &dsa, int nbc, const BCond *bcond,
	               int numCornerNodes, const std::vector<int> &cornerNodes, const std::vector<DofSet> &cornerDofs,
	               int ncbc = 0, const ComplexBCond *cbcond = nullptr, int numWetInterfaceNodes = 0,
	               const int *wetInterfaceNodes = nullptr, const DofSet *wetInterfaceDofs = nullptr);

	ConstrainedDSA(const DofSetArray &dsa, int n, const int *sing = nullptr);  // PJSA: 1-23-01
	ConstrainedDSA(const DofSetArray &dsa, const ConstrainedDSA &cdsa, int n, const int *sing = 0);

	ConstrainedDSA(const DofSetArray &dsa, const ConstrainedDSA &cdsa);

	ConstrainedDSA(ConstrainedDSA &&) = default;

	~ConstrainedDSA() override = default;

private:
	int invrowcolmax;
};

// Auxiliary definitions
#define BCFREE  0
#define BCLOAD  1
#define BCFIXED 2

class SimpleNumberer : public EqNumberer {
public:
	SimpleNumberer(int nnodes, int *renumb = nullptr, int myMap = 0);

	~SimpleNumberer() override = default;

	void setWeight(int, int);

	void makeOffset();
};

#endif
