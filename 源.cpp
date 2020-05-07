#include "mySturct.h"
#include"k-merCount.h"
#include <windows.h>
#include<time.h>
using namespace std;


void saveData()
{
    while (1)
    {
        if (isReadOver)
        {
            bool isCountOver = false;
            while (!isCountOver)
            {
                int i = 0;
                for (i; i < num_of_bitset; i++)
                {
                    if (myBitCode[i].isRead || !myBitCode[i].isWrite)
                        break;
                }
                if (i >= num_of_bitset)
                    isCountOver = true;
            }
            count_of_readHashSlotCall++;
            cout << "��ʼд����" << endl;

      
            string fileName = "temp" + to_string(count_of_readHashSlotCall) + ".txt";
            ofstream outfile(fileName, ios::app);
            //outfile << monitor_memUsage_call << endl;
            for (int i = 0; i < hash_size; i++)
            {
                readHashSlot(i, outfile);
            }
            cout << "�������н���" << endl;
            outfile.close();
            return;
        }
    }


}

int main()
{
    clock_t startTime, endTime;
    startTime = clock();//��ʱ��ʼ
    string readPath = "data\\hs_ref_GRCh38.p12_chr1.fa";
    thread readThread(readfile, readPath);
    readThread.detach();
    thread countThreads[num_of_count_thread];
    for (int i = 0; i < num_of_count_thread; i++)
        countThreads[i] = thread(countKmer);
    for (int i = 0; i < num_of_count_thread-1; i++)
        countThreads[i].detach();
    countThreads[num_of_count_thread - 1].join();

    bool isCountOver = false;
    while (!isCountOver)
    {
        int i = 0;
        for (i; i < num_of_bitset; i++)
        {
            if (myBitCode[i].isRead || !myBitCode[i].isWrite)
                break;
        }
        if (i >= num_of_bitset)
            isCountOver = true;
    }
    endTime = clock();//��ʱ����
    cout << "The count run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
   
    
    thread monitorMemUsage(saveData);
    monitorMemUsage.join();
    endTime = clock();//��ʱ����
    cout << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;


    
    system("pause");
    return 0;

}

