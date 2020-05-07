#ifndef K_MERCOUNT_H
#define K_MERCOUNT_H
#include "mySturct.h"
#include <string> 
#include<sstream>
#include <fstream> //文件操作库函数
#include <iostream> 
#include <thread>
#include <mutex>
using namespace std;



//A-00,C-01,G-10,T-11
const int temp_item_vector_size = 1024;
vector<read>myBitCode(num_of_bitset);
treeFactory instance;
bool isReadOver = false;
mutex isCountOverMut;
bool isCountOver = false;
bool hashAvil = true;
mutex hashAvilMut;
inline bool isDescriptionString(char t)
{
    if (t == '>')
        return true;
    return false;
}

//A-00,C-01,G-10,T-11
inline item countComplement(const item& it)
{
    item complement;
    int i = 0;
    for (int i = 0; i < k * size_of_codealph; i++)
    {
        if (it.kmer[i] == 0)
            complement.kmer[i] = 1;
        else
            complement.kmer[i] = 0;
    }
    return complement;
   
}
string binaryToAlph(string binaryStr)
{
    int i = 0;
    string res;
    while (i < binaryStr.size())
    {
        if (binaryStr[i] == '1')
        {
            i++;
            if (binaryStr[i] == '1')
            
                res += 'T';
            else
                res += 'G';
        }
        else
        {
            i++;
            if (binaryStr[i] == '1')
                res += 'C';
            else
                res += 'A';

        }
        i++;
    }
    return res;
}

void bitsetCpy(int start, int end, readUnit *source, kmerVector& des)//从readUnit中读取一个K-mer
{
    int count = 0;
    for (int i = start; i <= end; i++, count++)
        des[count] = source->bitCode[i];
}


void pushDataInHash(vector<item>& t)
{
    sort(t.begin(), t.end(), cmp);
    int i = 0;
    while (i < t.size())
    {
        int curIndex = i;
        int next = i + 1;
        string s1 = t[i].kmer.to_string();
        string s2;
        if (i < t.size() - 1)
            s2 = t[next].kmer.to_string();
        while (i < t.size() - 1&&s1==s2)
        {
            string s1 = t[i].kmer.to_string();
            ++t[curIndex];
            i++;
            next++;
        }
        int index_in_hash = t[curIndex].leftMove(hash_size_exp);
     //   string s = t[curIndex].kmer.to_string();
            
        if (res_Hash[index_in_hash] == NULL)
                res_Hash[index_in_hash] = new leafNode(-1, k * size_of_codealph - hash_size_exp, NULL);
        midNode* root = res_Hash[index_in_hash];
        leafNode* res_leaf = instance.searchNode(root, t[curIndex]);
        

        ostringstream oss;
        oss << std::this_thread::get_id();
        string stid = oss.str();
        long long tid = std::stoull(stid);
        long long expect_val = -1;

        while (!res_leaf->thread_index.compare_exchange_weak(expect_val, tid))
            continue;
       
        instance.addData(t[curIndex], res_leaf, index_in_hash);
        
        while (!res_leaf->thread_index.compare_exchange_weak(tid, expect_val))
            continue;


        i++;
    }

    //释放内存
    t.clear();
    t.shrink_to_fit();

}

int curnum_of_read = 0;
int curnum_of_line = 0;
mutex cur_countRead_mut;
void readfile(string readPath)
{
    //cout << "read" << endl;
    ifstream myfile(readPath);
    if (!myfile.is_open())
    {
        cout << "未成功打开文件" << endl;
        return ;
    }
   
    string temp;
    int index_of_bitset = 0;
    long long int maxCount = 0;
    long long int countOfsizeGreater1MB = 0;
    long long int countOfsizeGreater8MB = 0;
    long long int countOfsizeGreater16MB = 0;
    long long int countOfsizeGreater32MB = 0;
    long long int totalOfreads = 0;



    while (myfile.peek() != EOF)
    {
        while (isDescriptionString(myfile.peek()))//读取每个read的说明字符串
        {
            getline(myfile, temp);
            curnum_of_line++;
        }
       

        while (!myBitCode[index_of_bitset].isWrite)//选择一个可写的read
        {
            int t = index_of_bitset + 1;
            index_of_bitset = t %num_of_bitset;
        }
        myBitCode[index_of_bitset].isWrite = false;
        myBitCode[index_of_bitset].initRead();



        int count = 0;//计算当前read的大小
        while (myfile.peek() != EOF )//读取一个read
        {
            int index = 0;//指向下一个读取的字符在temp中的位置
            getline(myfile, temp);
            curnum_of_line++;
           // cout << temp << endl;

         

            while (index < temp.size())
            {

                if (myBitCode[index_of_bitset].cur->size <size_of_bitset)
                {
                    switch (temp[index])
                    {
                    case 'A':
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++, 0);//将myBitCode[index_of_bitset].cur->size位置上的值设为0
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++, 0);
                        count += 2;
                        break;
                        //A-00,C-01,T-11,G-10
                    case 'C':
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++, 0);
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++);
                        count += 2;
                        break;
                    case 'T':
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++);
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++);
                        count += 2;
                        break;
                    case 'G':
                        //A-00,C-01,G-10,T-11
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++);
                        myBitCode[index_of_bitset].cur->bitCode.set(myBitCode[index_of_bitset].cur->size++,0);
                        count += 2;
                        break;
                    default:
                        break;
                    }
                }
                else
                {
                    if (myBitCode[index_of_bitset].num_of_readUnit >= Maxnum_of_readUnit)
                    {
                        int next_index = 0;
                        while (!myBitCode[next_index].isWrite)//选择一个可写的read
                        {
                            int t = next_index + 1;
                            next_index = t % num_of_bitset;
                        }
                        myBitCode[next_index].isWrite = false;
                        myBitCode[next_index].initRead();
                        int t = (k - 1) * size_of_codealph;
                       // cout << myBitCode[index_of_bitset].cur->size << endl;;
                       // cout << index_of_bitset << endl;
                        int j = myBitCode[index_of_bitset].cur->size - t;
                        int i = 0;
                        for (i,j; i < t,j<myBitCode[index_of_bitset].cur->size; i++,j++)
                        {
                            myBitCode[next_index].cur->bitCode[i] = myBitCode[index_of_bitset].cur->bitCode[j];
                        }
                        myBitCode[next_index].cur->size = t;
                        if (i != t || j != myBitCode[index_of_bitset].cur->size)
                        {
                            cout << "更换read缓存出错" << endl;
                        }
                        curnum_of_read++;
                       // cout << "curnum_of_read:"<<curnum_of_read << endl;
                        myBitCode[index_of_bitset].isRead = true;
                        index_of_bitset = next_index;
                       // cout << index_of_bitset << endl;

                    }
                    else
                    {
                        myBitCode[index_of_bitset].addBitset();
                        myBitCode[index_of_bitset].num_of_readUnit++;
                    }
                    continue;
                }
                index++;
            }

            if (isDescriptionString(myfile.peek())|| myfile.peek() == EOF)//如果该read读取完毕
            {
                cur_countRead_mut.lock();
                curnum_of_read++;
                //cout << curnum_of_read << endl;
                cur_countRead_mut.unlock();
                myBitCode[index_of_bitset].isRead = true;
                break;
            }
            

        }

        //测试用
       // myBitCode[index_of_bitset].cur->size = 0;
        maxCount = count > maxCount ? count : maxCount;
        if (count > 1024 * 1024 * 8 * 32 - 1)
            countOfsizeGreater32MB++;
        else if (count > 1024 * 1024 * 8 * 16 - 1)
            countOfsizeGreater16MB++;
        else if (count > 1024 * 1024 * 8 * 8 - 1)
            countOfsizeGreater8MB++;
        else if (count > 1024 * 1024 * 8 - 1)
            countOfsizeGreater1MB++;
        totalOfreads++;

    }

    isReadOver = true;
    string writePath = "temp.txt";
    ofstream outfile(writePath, ios::app);

    outfile << "----------" << readPath << "的统计数据" << "----------" << endl;
    outfile << "maxCount  " << maxCount << endl << "countOfsizeGreater1MB  " << countOfsizeGreater1MB << endl;
    outfile << "countOfsizeGreater8MB  " << countOfsizeGreater8MB << endl;
    outfile << "countOfsizeGreater16MB  " << countOfsizeGreater16MB << endl;
    outfile << "countOfsizeGreater32MB  " << countOfsizeGreater32MB << endl;
    outfile << "totalOfreads  " << totalOfreads << endl;
    outfile.close();
    myfile.close();

}


int cur_countRead = 0;

void countKmer()
{
    while (1)
    {
        if (isCountOver)
            return;
        //cout << "count" << endl;
        vector<item> temp_item_vector;
        int index_of_count = 0;
        for (index_of_count; index_of_count < num_of_bitset; index_of_count++)
            if (myBitCode[index_of_count].isRead)
                break;
        if (index_of_count >= num_of_bitset)
        {
            if (!isReadOver)
                continue;
            else 
            {
                int i = 0;
                for(i;i<num_of_bitset;i++)
                    if (myBitCode[i].isRead)
                        break;
                if (i >= num_of_bitset)
                {
                    return;
                }
                     
            }
        }
        ostringstream oss;
        oss << std::this_thread::get_id();
        string stid = oss.str();
        long long tid = std::stoull(stid);
        long long expect_val = -1;

       
        if (myBitCode[index_of_count].thread_id.compare_exchange_weak(expect_val, tid))
        {
            myBitCode[index_of_count].isRead = false;
        }
        else
        {
            continue;
        }
       // cout << myBitCode[index_of_count].cur->size;
        //将当前read的cur指针设置为头指针
        myBitCode[index_of_count].cur = myBitCode[index_of_count].getHead();
        if (myBitCode[index_of_count].cur->size < k*2)//如果当前read的大小小于k-mer的大小，则直接返回
        {
            cout << myBitCode[index_of_count].cur->size;
            myBitCode[index_of_count].initRead();
            myBitCode[index_of_count].isRead = false;
            myBitCode[index_of_count].isWrite = true;
            continue;
        }

        int preindex = k*size_of_codealph - 1, rearindex = 0;
        readUnit* preBitset = myBitCode[index_of_count].cur;
        readUnit* rearBitset = myBitCode[index_of_count].cur;

        while (preBitset != NULL)
        {
            item temp;
            kmerVector tempKmer;
            item tempcomplement;
            string s;
            string sc;
            while (preindex < preBitset->size && rearindex < rearBitset->size)
            {
                bitsetCpy(rearindex, preindex, preBitset, tempKmer);
                temp.kmer = tempKmer;
                s = tempKmer.to_string();
                tempcomplement = countComplement(temp);
                sc = tempcomplement.kmer.to_string();
                if (s > sc)
                    temp_item_vector.push_back(temp);
                else
                    temp_item_vector.push_back(tempcomplement);


                if (temp_item_vector.size() > temp_item_vector_size)
                    pushDataInHash(temp_item_vector);

                preindex++;
                rearindex++;
            }
            if (preindex >= preBitset->size)
            {
                preindex = 0;
                preBitset = preBitset->next;
            }
            while (preBitset != NULL && preindex < preBitset->size && rearindex < rearBitset->size)
            {
                for (int i = rearindex; i < rearBitset->size; i++)
                    bitsetCpy(rearindex, rearBitset->size - 1, rearBitset, tempKmer);
                int preStart = rearBitset->size - rearindex + 1;
                for (int i = preStart, j = 0; i < k * size_of_codealph, j < preindex; i++, j++)
                {
                    tempKmer[i] = preBitset->bitCode[j];
                }

                temp.kmer = tempKmer;
                s = tempKmer.to_string();
                tempcomplement = countComplement(temp);
                sc = tempcomplement.kmer.to_string();
                if (s > sc)
                    temp_item_vector.push_back(temp);
                else
                    temp_item_vector.push_back(tempcomplement);


                if (temp_item_vector.size() > temp_item_vector_size)
                    pushDataInHash(temp_item_vector);

                preindex++;
                rearindex++;
            }
            if (rearindex >= rearBitset->size)
            {
                rearBitset = rearBitset->next;
                rearindex = 0;
            }
        }
        while (!myBitCode[index_of_count].thread_id.compare_exchange_weak(tid, expect_val))
        {
            continue;
        }
        cur_countRead_mut.lock();
        cur_countRead++;
        cout << "cur_countRead: " << curnum_of_read << endl;
        cur_countRead_mut.unlock();
        myBitCode[index_of_count].isWrite = true;
        if (temp_item_vector.size() > 0)
            pushDataInHash(temp_item_vector);
    }
    return;
}


static int count_of_readHashSlotCall = 0;

void saveleafNode(leafNode* Node, string& curMer, ofstream& outfile)
{
    Node->sortAndMergeDataUnit();
    int i = 0;
    while (i < Node->index)
    {
        int t = 0;
        string subMer;
        while (i < Node->index && t < Node->size)
        {
            if (Node->dataUnit[i])
                subMer += "1";
            else
                subMer += "0";
            t++;
            i++;
        }
        t = 0;
        bitset<size_of_itemCount> tempCount;
        while (i < Node->index && t < size_of_itemCount)
        {
            tempCount[t] = Node->dataUnit[i];
            t++;
            i++;
        }
        unsigned long count = tempCount.to_ulong();
        string Kmer = binaryToAlph(curMer + subMer);
        //cout << Kmer << "\t" << count << endl;
        outfile << Kmer << "\t" << count << endl;
  
    }
    return;
}
void readNode(midNode* Node,string &curMer, ofstream& outfile)
{
    if (Node == NULL)
        return;
    if (Node->id() == 0)
    {
        for (int i = 0; i < num_of_nodeChild; i++)
        {
            if (Node->childNode[i] == NULL)
                continue;
            switch (i)
            {
            case 0:
                curMer += "00";
                readNode(Node->childNode[i], curMer,outfile);
                break;
            case 1:
                curMer += "01";
                readNode(Node->childNode[i], curMer,outfile);
                break;
            case 2:
                curMer += "10";
                readNode(Node->childNode[i], curMer,outfile);
                break;
            case 3:
                curMer += "11";
                readNode(Node->childNode[i], curMer,outfile);
                break;
            default:
                break;
            }
           
        }
        if(curMer.size()>0)
            for (int i = size_nodeChild_code; i > 0; i--)
                curMer.pop_back();
        delete Node;
        Node = NULL;
        return;
    }
    else 
    {
        leafNode* leaf = reinterpret_cast<leafNode*>(Node);
        saveleafNode(leaf, curMer, outfile);
        delete Node;
        Node = NULL;
        for(int i= size_nodeChild_code;i>0;i--) 
            curMer.pop_back();
        return;
    }
}
void readHashSlot(int index, ofstream &outfile)
{
    if (res_Hash[index] == NULL)
        return;
    string curMer ;
    bitset<hash_size_exp> in;
    in.reset();
    midNode* root = res_Hash[index];
    int testIndex = index;
    int i = 0;
    while (index&&i<size_of_itemCount)
    {
        bool t = index % 2;
        if (t)
        {
            in[i] = 1;
        }
        i++;
        index=index >> 1;
    }
    curMer += in.to_string();
    //cout << testIndex << "\t" << curMer << endl;
    readNode(root, curMer, outfile);
    //delete res_Hash[index];
    return;
}

#endif // !K-MERCOUNT_H

