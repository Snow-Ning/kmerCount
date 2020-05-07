#ifndef MY_STURCT_H
#define MY_STURCT_H
#include <bitset>
#include <vector>
#include<mutex>
#include<iostream>
#include<algorithm>
#include<atomic>
#include <math.h>
using namespace std;

class leafNode;
class midNode;

const int k = 30;//k-mer��С

const int size_of_bitset = 1024*1024*8;//һ��readUnit�Ĵ�С
const short Maxnum_of_readUnit = 2;
const int num_of_bitset = 40;

const int num_of_count_thread = 32;

const int hash_size = 1024 * 1024;
const short hash_size_exp = 20;
const int size_of_codealph = 2;
const int size_of_itemCount = 20;

const int size_of_nodeUnit = 1024*(k*size_of_codealph-hash_size_exp+size_of_itemCount)/8;
const int num_of_nodeChild = 4;
const int size_nodeChild_code = 2;



midNode* res_Hash[hash_size];

typedef bitset<k*size_of_codealph> kmerVector;
typedef bitset<size_of_bitset> readUnitVector;
struct item;

vector<item> countGreaterMax;
mutex GreaterMaxMut;
struct item
{
public:
	kmerVector kmer;
	bitset<size_of_itemCount> count;

	item()
	{
		count[0] = 1;
	}
	item operator+(const item& t)
	{
		bool pro = 0;
		int i = 0;
		int num_of_1 = 0;
		while (i < size_of_itemCount)
		{
			num_of_1 = 0;
			if (this->count[i])
				num_of_1++;
			if (t.count[i])
				num_of_1++;
			if (pro)
				num_of_1++;

			switch (num_of_1)
			{
			case 0:
				this->count[i] = 0;
				pro = 0;
				break;
			case 1:
				this->count[i] = 1;
				pro = 0;
				break;
			case 2:
				this->count[i] = 0;
				pro = 1;
				break;
			case 3:
				this->count[i] = 1;
				pro = 1;
				break;
			default:
				cout << "There is a problem in the addition of item." << endl;
				break;
			
			}
			i++;
		}
		if (pro)
		{
			GreaterMaxMut.lock();
			cout << "There is a problem in the addition of item." << endl;
			countGreaterMax.push_back(t);
			for (int i = 0; i < size_of_itemCount; i++)
				count[i] = 0;
			GreaterMaxMut.unlock();
		}
		return *this;
	}
	item operator++()
	{
		int i = 0;
		for (i; i < size_of_itemCount; i++)
		{
			if (count[i])
				count[i] = 0;
			else if (!count[i])
			{
				count[i] = 1;
				return *this;
			}
		}
		if (i >= size_of_itemCount)
		{
			cout << "item ++ error" << endl;
			//system("pause");
		}
		return *this;
	}
	bool isequal(const item& t)
	{
		int i = 0;
		while (i < k * size_of_codealph)
		{
			if (this->kmer[i] != t.kmer[i])
				return false;
			i++;
		}
		return true;
	}
	int leftMove(int t)
	{
		int j = 0;
		int result = 0;
		string s = kmer.to_string();
		//cout << s << endl;
		int exp = 1;
		for (int i = 0; i < t; i++)
		{
			result += kmer[i]*exp;
			exp *= 2;
		}
		for (int i = t; i < kmer.size(); i++,j++)
		{
			kmer[j] = kmer[i];
		}
		for (int i = kmer.size() - t; i < kmer.size(); i++)
		{
			kmer[i] = 0;
		}
		return result;
	}
};

bool cmp(const item &x, const item &y)
{
	int t = 0;
	
	while (t < k * size_of_codealph)
	{
		
		if (x.kmer[t] != y.kmer[t])
		{
			if (x.kmer[t])
				return true;
			else
				return false;
		}
		t++;
	}
	return false;
}


struct readUnit
{
public:
	readUnitVector bitCode;
	int size=0;//bitCode��ʵ�ʰ����ı�����
	readUnit* next = NULL;
};

struct mytimer
{
	
};

class read 
{
private:
	readUnit* head = NULL;
public:
	short num_of_readUnit = 1;
	atomic<long long>thread_id = -1;
	readUnit * cur =NULL;
	bool isWrite = true;//���ڶ��߳�
	bool isRead = false;//���ڶ��߳�
	read()
	{
		head = NULL;
		cur = NULL;
		isWrite = true;
		isRead = false;
	}
	void addBitset()
	{
		readUnit* temp = new(readUnit);
		cur->next = temp;
		cur = cur->next;
	}
	readUnit* getHead()
	{
		return head;
	}
	void setHead()
	{
		head = new(readUnit);
	}
	void initRead()
	{
		if (head == NULL)
		{
			setHead();
		}
		while (head->next != NULL)
		{
			cur = head;
			head = cur->next;
			delete cur;
		}
		head->size = 0;
		num_of_readUnit = 1;
		cur = head;
		return;
	}
};

class midNode
{
public :
	midNode* childNode[num_of_nodeChild];
	virtual int id()
	{
		return 0;
	}
	virtual int getSize()
	{
		return 0;
	}
};


typedef bitset<size_of_nodeUnit> dataUnit_of_node;
class leafNode :public midNode
{
public :
	short pos_in_parent = 0;
	short height = 0;
	atomic<long long> thread_index = -1;
	int size = 0;//ָʾ��ǰ�߶ȵĽ��Ӧ�ô洢��λ�ĸ���
	int index = 0;//4 �ֽ�ָ��ýڵ��ڴ����һ���洢λ
	midNode* parent = NULL;//8Byte
	dataUnit_of_node  dataUnit;//128Byte

	leafNode()
	{

	}
	leafNode(const short &pos, const int &s, midNode* p)
	{
		pos_in_parent = pos;
		size = s;
		parent = p;
	}

	~leafNode()
	{
		~dataUnit_of_node();
		index = size = pos_in_parent = 0;
		parent = NULL;
	}

	virtual int id()
	{
		return 1;
	}
	virtual int getSize()
	{
		return size;
	}

	

	void sortAndMergeDataUnit()
	{
		vector<item> sortVector;
		int count = 0;
		while (count < index)
		{
			int t = 0;
			item tempItem;

			//Ϊitem�е�kmer��ֵ
			while (count < index && t <  size)
			{
				tempItem.kmer[t] = this->dataUnit[count];
				t++;
				count++;
			}
			//Ϊitem�е�count��ֵ
			t = 0;
			while (count < index && t < size_of_itemCount)
			{
				tempItem.count[t] = this->dataUnit[count];
				t++;
				count++;
			}
			sortVector.push_back(tempItem);
		}

		count = 0;
		//����ʱ�����Vector�е�item���򣬲��ϲ�kmer��ͬ����
		sort(sortVector.begin(), sortVector.end(), cmp);
		for (int i = 0; i < sortVector.size();)
		{
			int j = i + 1;
			while (j < sortVector.size() && sortVector[i].isequal(sortVector[j]))
			{
				unsigned long  x = sortVector[i].count.to_ulong();
				unsigned long y = sortVector[j].count.to_ulong();
				if (x + y >= pow(2, size_of_itemCount))
				{
					//cout <<"x: "<< sortVector[i].count.to_string() << "\t" << x << endl;
					//cout <<"y: " <<sortVector[j].count.to_string() << "\t" << x << endl;
					for (int m = 0; m < size_of_itemCount; m++)//��ǰһ��item���count��Ϊ���ֵ
						sortVector[i].count[m] = 1;
					x = x + y - pow(2, size_of_itemCount) + 1;
					sortVector[j].count=bitset<size_of_itemCount>(x);
					break;
				}
				else
				{
					//sortVector[i].count = bitset<size_of_itemCount>(x+y);
					sortVector[i] = sortVector[i] + sortVector[j];
				}
				j++;
			}

			//���Ѻϲ��õ�һ��kmerд�����ڴ�
			int t = 0;//��Ϊvector��item���±�
			//�洢item�е�kmer
			while (count < index && t < size )
			{
				this->dataUnit[count] = sortVector[i].kmer[t];
				t++;
				count++;
			}
			//�洢kmer�е�count
			t = 0;
			while (count < index && t < size_of_itemCount)
			{
				this->dataUnit[count] = sortVector[i].count[t];
				t++;
				count++;
			}
			i = j;
		}
		index = count;
		
	}

	//parent: �ýڵ�ĸ��ڵ�
	//pos: �ý���ڸ��ڵ����λ

	 bool push_back(item data)
	{
		int t = 0;
		while (index < size_of_nodeUnit&&t<size)//�洢k-mer
		{

			dataUnit[index] = data.kmer[t];
			t++;
			index++;	
		}
		t = 0;
		while (index < size_of_nodeUnit && t < size_of_itemCount)//�洢count
		{
			dataUnit[index] = data.count[t];
			t++;
			index++;
		}
		return true;
	};

	 
};

class secondLeafNode :public midNode
{
public:
	short pos_in_parent = 0;
	int index = 0;
	int size = 0;
	readUnit* head = NULL;
	readUnit* cur = NULL;
	void addBitset()
	{
		readUnit* temp = new(readUnit);
		cur->next = temp;
		cur = cur->next;
	}
	readUnit* getHead()
	{
		return head;
	}

};

class treeFactory
{
public:
	midNode* BurstNode(leafNode* leaf, const int& index_in_hash);
	leafNode* searchNode(midNode*, item&);
	bool  addData(item& it, leafNode* leaf,const int&index_in_hash)  //�ڽ���в���item
	{
		if (leaf == NULL)
			return false;
		int index = leaf->index;
		int size = leaf->size;
		if (index + size + size_of_itemCount >= size_of_nodeUnit)//�����ǰ�ڵ���ڴ治��
		{

			leaf->sortAndMergeDataUnit();
			index = leaf->index;
			if (index >= size_of_nodeUnit >> 1)
			{
				if(leaf->parent==NULL)

				midNode*temp=BurstNode(leaf,index_in_hash);
				

				return false;
			}
			leaf->push_back(it);
		}
		else
			leaf->push_back(it);
		return true;
			
	}


};


    // leaf Ҫ���ѵĽ��
	//1������һ���м���temp
	//2����leaf��dataUnit�����ݴ�ŵ�temp���ӽ����
	//3����leaf�ĸ��ڵ�ָ������ָ���Ϊָ��temp
	//4��delete leaf

	//1���ȼ����ӽ���ܷ��ٲ�������
	//2������ܣ�����leafNode��push_back�����������ݣ�������ܣ����Ѹý�㣬�ٲ�������
midNode* treeFactory::BurstNode(leafNode* leaf,const int &index_in_hash)//����һ��Ҷ�ӽ��
{
	midNode* temp = new midNode();//����һ���м���
	int count = 0;
	int size = leaf->size;
	int index = leaf->index;
	int childSize = size - size_nodeChild_code;

	//�����ӽڵ㲢Ϊ��dataUnit��ֵ
	while (count < index)
	{
		short t = 0;
		short pos = 0;
		while (count < index && t < size_nodeChild_code)//��ȡ�ò���itemӦ�÷����м�����ĸ��ֽڵ�
		{
			pos = pos * 2 + leaf->dataUnit[count];
			t++;
			count++;
		}
		if (temp->childNode[pos] == NULL)
		{
			temp->childNode[pos] = new leafNode(pos, childSize, temp);
		}

		//Ϊitem�е�kmer��count��ֵ
		t = 0;
		item tempItem;
		while (count < index && t < childSize)
		{
			tempItem.kmer[t] = leaf->dataUnit[count];
			t++;
			count++;
		}
		t = 0;
		while (count < index && t < size_of_itemCount)
		{
			tempItem.count[t] = leaf->dataUnit[count];
			t++;
			count++;
		}
		addData(tempItem, (leafNode*)temp->childNode[pos],index_in_hash);
	}
	if (leaf->parent == NULL)
	{
		res_Hash[index_in_hash] = temp;
		return temp;
	}
	leaf->parent->childNode[leaf->pos_in_parent] = temp;
	delete leaf;
	leaf = NULL;
	return temp;
}


int count_searchNode = 0;
leafNode* treeFactory::searchNode(midNode* Node, item& it)//����itemӦ�ò���Ľ��
{
	count_searchNode++;
	if (Node == NULL)
	{
		cout << count_searchNode << endl;
		cout << "there is nullptr in searchNode" << endl;
		system("pause");
		return NULL;
	}
	int id_Node = Node->id();
	if (id_Node == 0)
	{
		int t = it.leftMove(size_nodeChild_code);
		int size = 0;
		if (Node->childNode[t] == NULL)
		{
			int i = 0;
			for (i; i < num_of_nodeChild; i++)
			{
				if (Node->childNode[i] != NULL)
				{
					size = Node->childNode[i]->getSize();
					break;
				}
			}
			if (i >= num_of_nodeChild)
			{
				cout << "There are error in searchNode" << endl;
				system("pause");
				return NULL;
			}
			Node->childNode[t] = new leafNode(t,size-size_nodeChild_code,Node);
		}
		return searchNode(Node->childNode[t], it);
	}
	else if (id_Node == 1)
	{
		leafNode* leaf = reinterpret_cast<leafNode*>(Node);
		return leaf;
	}
}

#endif
