#ifndef _KDTREEH_
#define _KDTREEH_
#pragma once
#include<array>
#include<vector>
#include<unordered_set>
#include<utility>
#include<algorithm>
#include<queue>
#include<cmath>
#include<random>
using std::pair;
using std::vector;
using std::array;
using std::unordered_set;
using std::hash;
using std::priority_queue;
using std::fabs;
namespace K_DIMENSIONAL
{
	
	template <typename Element,unsigned Dimension, typename DimType=int>
	using Point = pair<Element, array<DimType, Dimension>>;
	template <typename Element, unsigned Dimension, typename DimType = int>
	using nearestNodeInfo = pair < double, Point<Element, Dimension, DimType>>;
	template <typename Element,unsigned Dimension, typename DimType=int>
	class KDTree
	{
	private:
		struct com
		{
			bool operator()(nearestNodeInfo<Element,Dimension,DimType> &n1, nearestNodeInfo<Element,Dimension,DimType> &n2)
			{
				return n1.first < n2.first;
			}
		};
	public:
		KDTree();
		KDTree(vector<Point<Element, Dimension, DimType>> &keyPoints);
		~KDTree();
		void clear();
		void createTree(vector<Point<Element, Dimension,DimType>> &keyPoints);
		bool insert(const Point<Element,Dimension,DimType> &p);
		bool remove(const Point<Element, Dimension, DimType> &p);
		priority_queue<nearestNodeInfo<Element, Dimension, DimType>, vector<nearestNodeInfo<Element, Dimension, DimType>>, com> findNearstNode(const array<DimType, Dimension> &point, unsigned number);
	private:
		struct kdNode
		{
			Point<Element, Dimension, DimType> dimData;
			int sortDim;
			kdNode *left;
			kdNode *right;
			kdNode *parent;
			kdNode() :left(nullptr), right(nullptr),parent(nullptr) {};
			kdNode(const Point<Element, Dimension, DimType> &data, int sd, kdNode *lt = nullptr, kdNode *rt = nullptr, kdNode *pt=nullptr)
				:dimData(data), sortDim(sd), left(lt), right(rt),parent(pt) {};
			kdNode &operator=(const kdNode &rhs)
			{
				dimData = rhs.dimData;
				sortDim = rhs.sortDim;
				left = rhs.left;
				right = rhs.right;
				parent = rhs.parent;
				return *this;
			}
		};
		struct priorityInfo
		{
			kdNode *ptr;
			double eucDistance;
			double verDistance;
			priorityInfo() :ptr(nullptr) {};
			priorityInfo(kdNode *p, double euc, double ver)
				:ptr(p), eucDistance(euc), verDistance(ver) {};
		};
		struct com1
		{
			bool operator()(const priorityInfo &p1, const priorityInfo &p2)
			{
				return p1.eucDistance > p2.eucDistance;
			}
		};
		struct mapping
		{
			size_t operator()(const Point<Element, Dimension, DimType> &p) const
			{
				size_t sum = 0;
				for (auto point : p.second)
					sum ^= hash<DimType>()(point);
				return hash<Element>()(p.first) ^ sum;
			}
		};
		kdNode *root;
		kdNode *nullNode;
		unsigned balanceCof;
		unordered_set<Point<Element, Dimension, DimType>,mapping> allPoints;
		void clear(kdNode *t);
		double calcDistance(const array<DimType, Dimension> &p1, const array<DimType, Dimension> &p2);
		unsigned findSortDim(const vector<Point<Element, Dimension, DimType>> &keyPoints);
		kdNode *findMidNode(vector<Point<Element, Dimension, DimType>> &keyPoints, int sortDim);
		kdNode *createTree(vector<Point<Element, Dimension, DimType>> &keyPoints, kdNode *pt);
		bool insert(const Point<Element, Dimension, DimType> &p, kdNode *&t,kdNode *pt,unsigned lastStep);
		bool remove(const array<DimType, Dimension> &deleteP, kdNode *&t);
	};
	template <typename Element, unsigned Dimension, typename DimType>
	KDTree<Element, Dimension, DimType>::KDTree<Element, Dimension, DimType>()
	{
		nullNode = new kdNode;
		nullNode->left = nullNode->right = nullNode;
		root = nullNode;
		vector<Point<Element, Dimension, DimType>> p;
		//以下为测试输入
		/*p.push_back({ 'A',{7,2} });
		p.push_back({ 'B',{5,4} });
		p.push_back({ 'C',{9,6} });
		p.push_back({ 'D',{2,3} });
		p.push_back({ 'E',{4,7} });
		p.push_back({ 'F',{8,1} });
		createTree(p);
		findNearstNode({2,4.5}, 2);*/
	}

	template<typename Element, unsigned Dimension, typename DimType>
	inline KDTree<Element, Dimension, DimType>::KDTree(vector<Point<Element, Dimension, DimType>>& keyPoints)
	{
		nullNode = new kdNode;
		nullNode->left = nullNode->right = nullNode;
		createTree(keyPoints);
	}

	template <typename Element, unsigned Dimension, typename DismType>
	KDTree<Element, Dimension, DismType>::~KDTree<Element, Dimension, DismType>()
	{
		clear();
		delete nullNode;
	}

	template<typename Element, unsigned Dimension, typename DimType>
	inline void KDTree<Element, Dimension, DimType>::clear()
	{
		clear(root);
		root = nullNode;
	}

	template<typename Element, unsigned Dimension, typename DimType>
	void K_DIMENSIONAL::KDTree<Element, Dimension, DimType>::createTree(vector<Point<Element, Dimension, DimType>>& keyPoints)
	{
		balanceCof = keyPoints.size() / 2;
		allPoints.clear();
		allPoints.insert(keyPoints.begin(), keyPoints.end());
		clear();
		root = createTree(keyPoints,nullptr);
	}

	template<typename Element, unsigned Dimension, typename DimType>
	bool K_DIMENSIONAL::KDTree<Element, Dimension, DimType>::insert(const Point<Element, Dimension, DimType>& p)
	{
		if (allPoints.find(p) == allPoints.end())
		{
			if (--balanceCof == 0)
			{
				allPoints.insert(p);
				balanceCof = allPoints.size() / 2;
				vector<Point<Element, Dimension, DimType>> temp(allPoints.begin(), allPoints.end());
				root = createTree(temp,nullptr);
				return true;
			}
			else
				return insert(p, root, nullptr, 0);
		}
		else
			return false;
		
	}

	template<typename Element, unsigned Dimension, typename DimType>
	bool K_DIMENSIONAL::KDTree<Element, Dimension, DimType>::remove(const Point<Element, Dimension, DimType> &p)
	{
		if (allPoints.find(p) != allPoints.end())
		{
			if (--balanceCof == 0)
			{
				allPoints.erase(p);
				balanceCof = allPoints.size() / 2;
				vector<Point<Element, Dimension, DimType>> temp(allPoints.begin(), allPoints.end());
				root = createTree(temp, nullptr);
				return true;
			}
			else
				return remove(p.second, root);
		}
		else
			return false;
	}

	/*                                      B     B     F       算     法                     */
	template<typename Element, unsigned Dimension, typename DimType>
	inline std::priority_queue<nearestNodeInfo<Element, Dimension, DimType>, vector<nearestNodeInfo<Element, Dimension, DimType>>,typename KDTree<Element,Dimension,DimType>::com> KDTree<Element, Dimension, DimType>::findNearstNode(const array<DimType, Dimension>& point, unsigned number)
	{
		/*Queue队列保存离确定点最近的number个点，按照距离最远的优先级最高的的方式构造该队列*/
		priority_queue<nearestNodeInfo<Element, Dimension, DimType>, vector<nearestNodeInfo<Element, Dimension, DimType>>, com> Queue;
		Queue.push({ FLT_MAX,Point<Element,Dimension,DimType>() });   //首先插入一个逻辑上距离point最远的空结点
		if (root == nullNode)
			return priority_queue<nearestNodeInfo<Element, Dimension, DimType>, vector<nearestNodeInfo<Element, Dimension, DimType>>, typename KDTree<Element, Dimension, DimType>::com>();
		kdNode *p = root;
		/*priQueu队列保存可能的回溯结点，按照距离最近的优先级最高的方式构造该队列*/
		priority_queue<priorityInfo, vector<priorityInfo>, com1> priQueue;
		priQueue.push(priorityInfo(p, calcDistance(point, p->dimData.second), fabs(point[p->sortDim] - p->dimData.second[p->sortDim])));  //首先将root结点插入当中
		int t = 0;   //逻辑迭代时间次数，保证在一定迭代次数后BBF算法可以退出
		while (!priQueue.empty())
		{
			t++;
			priorityInfo temp = priQueue.top();
			priQueue.pop();
			/*如果Queue的最高优先级项的邻近距离小于查找点到当前点确定的分隔超平面的距离则不访问该点的分支*/
			if (Queue.top().first < temp.verDistance)
				continue;
			kdNode *q = temp.ptr;
			while (q!=nullNode)
			{
				t++;
				/*距离当前结点的欧式距离*/
				double dis = calcDistance(point, q->dimData.second);
				/*判断是否需要更新k最邻近距离*/
				if (dis < Queue.top().first)
				{
					if (Queue.size() < number)
						Queue.push({ dis,q->dimData });
					else
					{
						Queue.pop();
						Queue.push({ dis,q->dimData });
					}
				}
				int sd = q->sortDim;
				/*循环访问后续结点，并将访问路径上的兄弟结点添加到PrQueue队列中*/
				if (point[sd] <= q->dimData.second[sd])
				{
					if (q->left != nullNode)
					{
						if (q->right != nullNode)
						{
							double distance = calcDistance(point, q->right->dimData.second);
							int st = q->right->sortDim;
							priQueue.push(priorityInfo(q->right, distance, fabs(point[st] - q->right->dimData.second[st])));
						}
						q = q->left;
					}
					else
						break;
				}
				else
				{
					if (q->right != nullNode)
					{
						if (q->left != nullNode)
						{
							double distance = calcDistance(point, q->left->dimData.second);
							int st = q->left->sortDim;
							priQueue.push(priorityInfo(q->left, distance, fabs(point[st] - q->right->dimData.second[st])));
						}
						q = q->right;
					}
					else
						break;
				}
			}
			if (t > 300)  //保证在300次迭代后BBF算法必须退出，即是可能未找到最优结点也要退出
				break;
		}
		return Queue;
	}

	template<typename Element, unsigned Dimension, typename DimType>
	inline void KDTree<Element, Dimension, DimType>::clear(kdNode * t)
	{
		if (t == nullNode);
		else
		{
			clear(t->left);
			clear(t->right);
			delete t;
		}
	}

	template<typename Element, unsigned Dimension, typename DimType>
	inline double KDTree<Element, Dimension, DimType>::calcDistance(const array<DimType, Dimension>& p1, const array<DimType, Dimension>& p2)
	{
		double distance = 0;
		if (p1.size() != p2.size())
			throw std::length_error("两点不匹配");
		else
		{
			for (unsigned i = 0; i < p1.size(); i++)
				distance += pow(p1[i] - p2[i], 2);
			return sqrt(distance);
		}
	}

	template<typename Element, unsigned Dimension, typename DimType>
	unsigned K_DIMENSIONAL::KDTree<Element, Dimension, DimType>::findSortDim(const vector<Point<Element, Dimension, DimType>>& keyPoints)
	{
		unsigned sortDim = 0;
		double sortVar = 0;
		for (unsigned i = 0; i < Dimension; i++)
		{
			double sum = 0;
			double avergae=0;
			double variance=0;
			for(auto point:keyPoints)
				sum += point.second[i];
			avergae = sum / keyPoints.size();
			for(auto point:keyPoints)
				variance += pow(point.second[i] - avergae, 2);
			variance /= keyPoints.size();
			if (sortVar < variance)
			{
				sortVar = variance;
				sortDim = i;
			}
		}
		return sortDim;
	}

	template<typename Element, unsigned Dimension, typename DimType>
	typename KDTree<Element, Dimension, DimType>::kdNode * K_DIMENSIONAL::KDTree<Element, Dimension, DimType>::findMidNode(vector<Point<Element, Dimension, DimType>>& keyPoints, int sortDim)
	{
		int mid = keyPoints.size() / 2;
		sort(keyPoints.begin(), keyPoints.end(), [sortDim](const Point<Element,Dimension,DimType> &p1,const Point<Element, Dimension, DimType> &p2) {return p1.second[sortDim] < p2.second[sortDim]; });
		kdNode *temp = new kdNode(keyPoints[mid], sortDim,nullNode,nullNode);
		return temp;
	}
	template<typename Element, unsigned Dimension, typename DimType>
	inline typename KDTree<Element, Dimension, DimType>:: kdNode * KDTree<Element, Dimension, DimType>::createTree(vector<Point<Element, Dimension, DimType>>& keyPoints, kdNode *pt)
	{
		if (keyPoints.empty())
			return nullNode;
		int sortDim = findSortDim(keyPoints);
		kdNode *temp = findMidNode(keyPoints, sortDim);
		vector<Point<Element, Dimension, DimType>> leftKeyPoints(keyPoints.begin(), keyPoints.begin() + keyPoints.size() / 2);
		vector<Point<Element, Dimension, DimType>> rightKeyPoints(keyPoints.begin() + keyPoints.size() / 2 + 1, keyPoints.end());
		temp->parent = pt;
		temp->left = createTree(leftKeyPoints,temp);
		temp->right = createTree(rightKeyPoints,temp);
		return temp;
	}
	template<typename Element, unsigned Dimension, typename DimType>
	inline bool KDTree<Element, Dimension, DimType>::insert(const Point<Element, Dimension, DimType>& p, kdNode *& t,kdNode *pt,unsigned lastStep)
	{
		static bool flag = false;
		if (t == root)
			flag = false;
		if (t == nullNode)
		{
			flag = true;
			int sortDim = static_cast<int>(lastStep + 1 - Dimension) >= 0 ? lastStep + 1 - Dimension : lastStep + 1;
			t = new kdNode(p, sortDim, nullNode, nullNode, pt);
		}
		else
		{
			if (p.second[t->sortDim] <= t->dimData.second[t->sortDim])
				insert(p, t->left, t, t->sortDim);
			else if (p.second[t->sortDim] > t->dimData.second[t->sortDim])
				insert(p, t->right, t, t->sortDim);
		}
		if (t == root)
			return flag;
	}
	template<typename Element, unsigned Dimension, typename DimType>
	inline bool KDTree<Element, Dimension, DimType>::remove(const array<DimType, Dimension>& deleteP, kdNode *&t)
	{
		static bool falg = false;
		static std::default_random_engine e;
		static std::bernoulli_distribution u;
		if (t == root)
			falg = false;
		if (t == nullNode)
			return false;
		else if (deleteP[t->sortDim] <= t->dimData.second[t->sortDim])
		{
			if(deleteP!=t->dimData.second)
				remove(deleteP, t->left);
			else if (t->left != nullNode && t->right != nullNode)
			{
				if (u(e))
				{
					t->dimData = t->left->dimData;
					remove(t->dimData.second, t->left);
				}
				else
				{
					t->dimData = t->right->dimData;
					remove(t->dimData.second, t->right);
				}
			}
			else
			{
				falg = true;
				kdNode *temp = t;
				t = t->right == nullNode ? t->left : t->right;
				delete temp;
			}
		}
		else if (deleteP[t->sortDim] > t->dimData.second[t->sortDim])
		{
			remove(deleteP, t->right);
		}
		if (t == root)
			return falg;
	}
}
#endif

