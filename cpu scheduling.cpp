
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <queue>
#include <functional>
#include <iomanip>
#include <thread>
#include <cmath>

using namespace std;


struct Process {
    char name;
    int arrival_time;
    int service_time;
    int remaining_time;
    int start_time;
    int finish_time;
    int turnaround;
    bool isCompleted;
    vector<int> isRunning;
    int state_size; 
    int priority;
    bool justfinished;
    
    

// Default constructor
    Process()
       : name('-'), arrival_time(0), service_time(0), remaining_time(0), start_time(-1), finish_time(-1),isCompleted(0),priority(0),justfinished(0) {}


    Process(char n, int arr_time, int srv_time,int state_size)
        : name(n), arrival_time(arr_time), service_time(srv_time), remaining_time(srv_time), start_time(-1), finish_time(-1),turnaround(-1),isRunning(state_size, -1) ,priority(0){}
            //state 0 --> process waiting
            //state 1--> process running
            //state -1 ---> process is idle (either finished or just created)   
};

string visualization , line;
int last_instant , n_of_processes,current_time=0;   


struct priorityComparator {
    bool operator()(const Process* p1, const Process* p2) {
        if( p1->priority > p2->priority )
            return true;
        else if((p1->priority == p2->priority) && ((p1->priority - p1->service_time) > (p2->priority - p2->service_time)))
            return true ;
        return false;

    }
};
struct ArrivalTimeComparator {
    bool operator()(const Process& p1, const Process& p2) {
        return p1.arrival_time < p2.arrival_time;  // Process with earlier arrival time has higher priority
    }
};
struct SPNComparator {
    bool operator()(const Process* p1, const Process* p2) {
        return p1->service_time < p2->service_time;  // Process with less service time has higher priority
    }
};

struct SRTComparator{
    bool operator()(const Process* p1, const Process* p2) {
        return p1->remaining_time < p2->remaining_time;  // Process with less service time has higher priority
    }
};

/*void moveHighestPriorityToFront(vector<Process*>& ready_queue) {
    if (ready_queue.empty()) {
        return; // Do nothing if the ready queue is empty
    }

    // Find the process with the highest priority
    int highest_priority_index = 0;
    for (size_t i = 1; i < ready_queue.size(); ++i) {
        if (ready_queue[i]->priority > ready_queue[highest_priority_index]->priority) {
            highest_priority_index = i;
        }
    }


    // Swap the process with the highest priority to the front
    if (highest_priority_index != 0) {
        swap(ready_queue[0], ready_queue[highest_priority_index]);
    }
}*/
void printCentered(const string& str,int isInt) {
    int len = str.length();
    int leftPadding;
    if(isInt)
        leftPadding = floor((5 - len) / 2.0);
    else
        leftPadding = ceil((5 - len) / 2.0);
    int rightPadding = 5 - len - leftPadding;
    cout << setw(leftPadding) << "" << str << setw(rightPadding) << "";
}

void printCenteredInt(int value) {
    string str = to_string(value);
    printCentered(str,1);
}

void printCenteredFloat(float value, int precision = 2) {
    ostringstream oss;
    oss << fixed << setprecision(precision) << value;
    printCentered(oss.str(),0);
}

bool inQueue(Process& p ,vector<Process*> &ready_queue){
    if (ready_queue.empty()) //if no queues were initialized yet
        return false;
    for (const Process *ready_p : ready_queue) {
        if (ready_p->name == p.name) {
            return true;
        }
    }
    return false;
}

bool multiQueueSearch(Process &p,vector<vector<Process*>> &ready_queues){
    if (ready_queues.empty()) //if no queues were initialized yet
        return false;

    for(vector<Process*> &ready_queue : ready_queues){
        if(inQueue(p,ready_queue))
            return true;
    }
    return false;

}

void updateToWaitQueue(vector<Process*>&ready_queue, Process &p,int time){
    if (!ready_queue.empty()) {
        for (Process *ready_p : ready_queue) {
                if (ready_p->name != p.name && !ready_p->isCompleted){
                    ready_p->isRunning[time]=0;
                    //cout<<"process "<<ready_p->name<<" is waiting at time "<< time<<endl;
                }
            } 
    }
}

void updateToWait(vector<vector<Process*>> &ready_queues, Process &p,int time){
    if(!ready_queues.empty()){
        for(vector ready_queue: ready_queues){
            updateToWaitQueue(ready_queue,p,time);
        }
    }
}

bool isNextQueuesEmpty(int queue_index,vector<vector<Process*>> &ready_queues){
    for(size_t i=queue_index+1;i<ready_queues.size();i++){
        if(!ready_queues[i].empty()){
            return false;
        }

    }
    return true;
}
bool isPrecedentQueueEmpty(int queue_index,vector<vector<Process*>> &ready_queues){
    for(int i=0;i<queue_index;i++){
        if(!ready_queues[i].empty()){
            return false;
        }

    }
    return true;

}
bool isReadyQueueSortedByArrivalTime(const vector<Process*>& ready_queue) {
    // Iterate through the ready queue and check if each process' arrival time is greater than or equal to the previous one
    for (size_t i = 1; i < ready_queue.size(); ++i) {
        if (ready_queue[i]->arrival_time < ready_queue[i - 1]->arrival_time) {
            return false; // If any process' arrival time is less than the previous, return false
        }
    }
    return true; // If all processes are in non-decreasing order, return true
}


void trace(vector<Process>& processes,int total_time, string scheduling_algorithm) {
    cout << scheduling_algorithm <<setw(6-scheduling_algorithm.size()+1);

    //print time
    for (int t = 0; t <= total_time; ++t) {
        cout << t % 10 << " ";
    }
    cout << endl;
    cout << "------------------------------------------------" << endl;


    for (int i = 0; i < n_of_processes; i++) {
        cout << processes[i].name<<setw(6)<< "|";
        for (int t = 0; t < total_time; t++) {
            if(t >= processes[i].arrival_time && t<processes[i].finish_time){
                if (processes[i].isRunning[t]==1) 
                    cout << "*|";  // Process is running
                else 
                    cout << ".|";  // Process is waiting  
            }
            else 
                cout << " |";  // Process is idle     
        }
        cout <<" "<< endl;
    }

    cout << "------------------------------------------------" << endl;
    cout<<endl;
    
}


void stats(vector<Process>& processes,string scheduling_algorithm){
    
    cout<<scheduling_algorithm<< endl;

    //display names of processes
    cout<< "Process" <<"    |";
    for (int i=0;i<n_of_processes;i++){
        cout<<"  "<<processes[i].name<<"  |";
    }
    cout<<endl;

    //display arrival time of processes
    cout<< "Arrival" <<"    |";
    for (int i=0;i<n_of_processes;i++){
        printCenteredInt(processes[i].arrival_time);
        cout<<"|";
    }
    cout<<endl;

    //display service time of processes
    cout<< "Service" <<"    |";
    for (int i=0;i<n_of_processes;i++){
        printCenteredInt(processes[i].service_time);
        cout<<"|";
    }
    cout<<" Mean|"<<endl;

    //display finish time of processes
    cout<< "Finish" <<"     |";
    for (int i=0;i<n_of_processes;i++){
       printCenteredInt(processes[i].finish_time);
        cout<<"|";
    }
    cout<<"-----|"<<endl;

    //calculate turnaround & avg_turnaround
    float avg_turnaround=0;
    cout<< "Turnaround" <<" |";
    for (int i=0;i<n_of_processes;i++){
       printCenteredInt(processes[i].turnaround);
        cout<<"|";
        avg_turnaround+=processes[i].turnaround;
    }
    avg_turnaround/=n_of_processes;
    printCenteredFloat(avg_turnaround);
    cout<<"|";
    cout<<endl;

    

    //calculate normturnaround & avg_normturnaround
    float avg_norm_turnaround=0,norm_turnaround,sum_norm_turnaround=0;
    cout<< "NormTurn" <<"   |";
    for (int i=0;i<n_of_processes;i++){
        norm_turnaround = processes[i].turnaround*1.0/processes[i].service_time;
        printCenteredFloat(norm_turnaround);
        cout<<"|";
        sum_norm_turnaround+=norm_turnaround;
    }

    //calculate the avg of normturnaround
    avg_norm_turnaround=sum_norm_turnaround/n_of_processes;
    printCenteredFloat(avg_norm_turnaround);
    cout<<"|";
    cout<<endl;
    cout<<endl;

}


void SPN (vector<Process> processes,int last_instant,bool trace_visualization){
    int i;
    int current_time = 0;
    vector<Process *> ready_queue;

    sort(processes.begin(), processes.end(), ArrivalTimeComparator());//sort processses by arrival time
    
    // Process each process in the order of shortest service time
    while(current_time < last_instant){

        //add processses that have arrived by the current time and isnt completed
        for(int j=0;j<n_of_processes;j++){
            Process *p = &processes[j];

            if (p->arrival_time <= current_time && !p->isCompleted && !inQueue(*p,ready_queue)){
                    ready_queue.push_back(p);
                    //update state to waiting state
                    p->isRunning[current_time]=0;
            }

        }
        // if queue is empty and there no processes have arrived yet, go to next instant
        if(ready_queue.empty())
            current_time++;


        //sort processes in the ready queue by shortest service time
        sort(ready_queue.begin(), ready_queue.end(), SPNComparator());

        if (!ready_queue.empty()){
            //get process with shortest service time & remove it from ready queue
            Process *selected_process = ready_queue.front();
            ready_queue.erase(ready_queue.begin());

            //update starting time of the process
            selected_process->start_time = current_time;

            //update current_time
            current_time += selected_process->service_time;

            //calculate time data for the process
            selected_process->finish_time = current_time;
            selected_process->turnaround = selected_process->finish_time - selected_process->arrival_time;

            //update state array of the process from arrival time to finish time
            for(i=selected_process->arrival_time;i<selected_process->finish_time;i++){

                if(i >= selected_process->start_time)
                    selected_process->isRunning[i]=1;//running state
                else 
                    selected_process->isRunning[i]=0;//waiting state
                
            }   
            //mark process as completed
            selected_process->isCompleted = true;
        }
       

    }


    if (trace_visualization) 
        trace(processes, last_instant, "SPN");
    else
        stats(processes,"SPN");

}

void FB1(vector<Process> processes,int last_instant,bool trace_visualization){
    vector<vector<Process*>> ready_queues; // Vector of dynamically allocated ReadyQueue pointers
    int number_of_queues =1,new_process=0;
    int i=0,flag ,current_time = 0;
    Process *selected_process;
    char last_running_process;
    
    ready_queues.push_back(vector<Process*>());//first ready queue with highest priority
    sort(processes.begin(), processes.end(), ArrivalTimeComparator());//sort processses by arrival time

    //iterate over the whole time interval
    while(current_time<last_instant)
    {
        //no process were chosen to run
        flag =0;
       //insert all processes that have arrived in the first ready queue
        for(int j=0;j<n_of_processes;j++){
            Process *p = &processes[j];
            if (p->arrival_time <= current_time && !p->isCompleted && !multiQueueSearch(*p,ready_queues)){
                ready_queues[0].push_back(p);
            }
        }
       
        while(i<number_of_queues ){
            //if the next queue contains the same process that started running
            if(!ready_queues[i].empty() && last_running_process == ready_queues[i].front()->name && !isNextQueuesEmpty(i,ready_queues))
                    i++;
            //if new process was admitted start from the first queue
            if(!ready_queues[0].empty()){
                i=0;
            }
            if(!ready_queues[i].empty()){
                //get first process in queue
                selected_process = ready_queues[i].front();
                ready_queues[i].erase(ready_queues[i].begin());

                //if it's its first time running, update start time
                if(selected_process->start_time == -1)
                    selected_process->start_time=current_time;

                //update state array & remaining time
                selected_process->isRunning[current_time]=1;
                selected_process->remaining_time-=1;

                //update the rest of processes to waiting state
                updateToWait(ready_queues,*selected_process,current_time);
                
                //if process hasnt completed task, insert in the next queue
                if(selected_process->remaining_time!=0)
                {
                    //admit new processes coming 
                    for(int j=0;j<n_of_processes;j++){
                        Process *p = &processes[j];
                        if (p->arrival_time == current_time+1 && !p->isCompleted && !multiQueueSearch(*p,ready_queues)){
                            new_process=1;
                        }
                    }
                    //if after admitting there is no change and all queues are empty, re-enter same queue
                    if(isPrecedentQueueEmpty(i,ready_queues) && isNextQueuesEmpty(i,ready_queues) && ready_queues[i].empty() && !new_process){
                        ready_queues[i].push_back(selected_process);
                    }
                    else{ 
                        if(i+1 < number_of_queues ){
                            ready_queues[i+1].push_back(selected_process);
                            
                        }
                        //create new queue if it's the last queue
                        else{
                            vector<Process*> new_queue = vector<Process*>();//create new queue
                            new_queue.push_back(selected_process); //insert process to new queue with lower priority
                            ready_queues.push_back(new_queue);//insert the newly created queue into the vector of queues
                            number_of_queues++;         
                        }
                    }
                }
                else
                {
                    selected_process->isCompleted=1;
                    selected_process->finish_time = current_time+1;
                    selected_process->turnaround = selected_process->finish_time - selected_process->arrival_time;
                }
                //mark flag indicating a process was run
                flag =1;
            }
            //only go to the next queue if this one is empty and there are more queues
            if(i < number_of_queues-1 && ready_queues[i].empty()){
                i++;
            }
            // if all queues are empty and there no processes have arrived yet, go to next instant
            if(isNextQueuesEmpty(i,ready_queues) && isPrecedentQueueEmpty(i,ready_queues) && ready_queues[i].empty())
                current_time++;

            //if a process was running, break to next instant
            if(flag){
                last_running_process = selected_process->name;
                break;
            }

        }
        new_process=0;
        current_time++;

    }
    
    if (trace_visualization) 
        trace(processes, last_instant, "FB-1");
    else
        stats(processes,"FB-1");

}

void FB2i(vector<Process> processes,int last_instant,bool trace_visualization){
    vector<vector<Process*>> ready_queues; // Vector of dynamically allocated ReadyQueue pointers
    int number_of_queues =1;
    int i=0,flag ,current_time = 0,start_time=0,new_process=0;
    Process *selected_process;
    char last_running_process;
    
    ready_queues.push_back(vector<Process*>());//first ready queue with highest priority
    sort(processes.begin(), processes.end(), ArrivalTimeComparator());//sort processses by arrival time

    //iterate over the whole time interval
    while(current_time<last_instant)
    {
        //no process were chosen to run
        flag =0;
       //insert all processes that have arrived in the first ready queue
        for(int j=0;j<n_of_processes;j++){
            Process *p = &processes[j];
            if (p->arrival_time <= current_time && !p->isCompleted && !multiQueueSearch(*p,ready_queues)){
                ready_queues[0].push_back(p);
            }
        }
       
        while(i<number_of_queues ){
            //if the next queue contains the same process that started running
            if(!ready_queues[i].empty() && last_running_process == ready_queues[i].front()->name && !isNextQueuesEmpty(i,ready_queues))
                i++;
            //if new process was admitted start from the first queue
            if(!ready_queues[0].empty()){
                i=0;
            }
            if(!ready_queues[i].empty()){
                //get first process in queue
                selected_process = ready_queues[i].front();
                ready_queues[i].erase(ready_queues[i].begin());

                //if it's its first time running, update start time
                if(selected_process->start_time == -1)
                    selected_process->start_time=current_time;

                //run process with specified quantum
                start_time=current_time;
                while(selected_process->remaining_time!=0 && (current_time - start_time < pow(2,i))){                    
                    //mark flag indicating a process was run
                    flag =1;
                    //update state array & remaining time
                    selected_process->isRunning[current_time]=1;
                    selected_process->remaining_time-=1;
                    //update the rest of processes to waiting state
                    updateToWait(ready_queues,*selected_process,current_time);
                    current_time++;
                }
               
                

                //if process hasnt completed task, insert in the next queue
                if(selected_process->remaining_time!=0)
                {
                    //admit new processes comming 
                    for(int j=0;j<n_of_processes;j++){
                        Process *p = &processes[j];
                        if (p->arrival_time == current_time && !p->isCompleted && !multiQueueSearch(*p,ready_queues)){
                            new_process=1;
                        }
                    }
                    //if after admitting there is no change and all queues are empty, re-enter same queue
                    if(isPrecedentQueueEmpty(i,ready_queues) && isNextQueuesEmpty(i,ready_queues) && ready_queues[i].empty() && !new_process){
                        ready_queues[i].push_back(selected_process);
                    }
                    else{
                        if(i+1 < number_of_queues ){
                            ready_queues[i+1].push_back(selected_process);

                        }
                        //create new queue if it's the last queue
                        else{
                            vector<Process*> new_queue = vector<Process*>();//create new queue
                            new_queue.push_back(selected_process); //insert process to new queue with lower priority
                            ready_queues.push_back(new_queue);//insert the newly created queue into the vector of queues
                            number_of_queues++;
                                    
                        }
                    }
                }
                else
                {
                    selected_process->isCompleted=1;
                    selected_process->finish_time = current_time;
                    selected_process->turnaround = selected_process->finish_time - selected_process->arrival_time;
                }
                
            }
            //only go to the next queue if this one is empty and there are more queues
            if(i < number_of_queues-1 && ready_queues[i].empty()){
                i++;
            }
            // if all queues are empty and there no processes have arrived yet, go to next instant
            if(isNextQueuesEmpty(i,ready_queues) && isPrecedentQueueEmpty(i,ready_queues) && ready_queues[i].empty())
                current_time++;

            //if a process was running, break to next instant
            if(flag){
                last_running_process = selected_process->name;
                break;
            }
            new_process=0;
        }

    }
    
    if (trace_visualization) 
        trace(processes, last_instant, "FB-2i");
    else
        stats(processes,"FB-2i");

}

void SRT(vector<Process> processes,int last_instant,bool trace_visualization){
    vector<Process*> ready_queue; // Vector of dynamically allocated ReadyQueue pointers
    int current_time=0;
    sort(processes.begin(), processes.end(), ArrivalTimeComparator());//sort processses by arrival time

    while (current_time < last_instant){
        //admit new processes
        for(int j=0;j<n_of_processes;j++){
                Process *p = &processes[j];
                if (p->arrival_time <= current_time && !p->isCompleted && !inQueue(*p,ready_queue)){
                    ready_queue.push_back(p);
                }
            }
       
        if(!ready_queue.empty() ){
            //sort by shortest remaining time
            sort(ready_queue.begin(), ready_queue.end(), SRTComparator());

            Process *selected_process = ready_queue.front();

        
            //remove the process
            ready_queue.erase(ready_queue.begin());
            //update start time if it's the first time running
            if(selected_process->start_time == -1)
                selected_process->start_time = current_time;

            updateToWaitQueue(ready_queue,*selected_process,current_time);

            //update state
            selected_process->isRunning[current_time]=1;
            selected_process->remaining_time -=1;

            //if process finished task update state & do not re-admit
            if(selected_process->remaining_time == 0){
                selected_process->finish_time=current_time+1;
                selected_process->isCompleted=1;
                selected_process->turnaround = selected_process->finish_time - selected_process->arrival_time;

            }
            else{
                ready_queue.push_back(selected_process);
            }
        }
        current_time++;

    }

     if (trace_visualization) 
        trace(processes, last_instant, "SRT");
    else
        stats(processes,"SRT");

}

void FCFS(vector<Process> processes, int last_instant, bool trace_visualization) {
    //  Sort processes by arrival time
    sort(processes.begin(), processes.end(), ArrivalTimeComparator());

    int current_time = 0;

    //  execute each process in arrival order
    for (auto& p : processes) {
        // wait if the current time is earlier than the process arrival
        if (current_time < p.arrival_time) {
            current_time = p.arrival_time;
        }

        // update time and process properties
        p.start_time = current_time;
        p.finish_time = current_time + p.service_time;
        p.turnaround = p.finish_time - p.arrival_time;

        // update isRunning state
        for (int t = p.arrival_time; t < p.finish_time; ++t) {
            if (t >= p.start_time) {
                p.isRunning[t] = 1; // Running
            } else {
                p.isRunning[t] = 0; // Waiting
            }
        }

        // move current time forward
        current_time = p.finish_time;

        // mark as completed
        p.isCompleted = true;
    }

    
    if (trace_visualization) {
        trace(processes, last_instant, "FCFS");
    } else {
        stats(processes, "FCFS");
    }
}
void RoundRobin(vector<Process> processes, int last_instant, int time_quantum, bool trace_visualization) {
    sort(processes.begin(), processes.end(), ArrivalTimeComparator()); // Sort by arrival time
    vector<Process*> ready_queue;
    int current_time = 0;
    int completed_processes = 0;
    
    

    while (completed_processes < n_of_processes || !ready_queue.empty()) {
    // Add processes that have arrived by the current time to the ready queue
        for (Process& p : processes) {
            if (p.arrival_time <= current_time && !p.isCompleted && !inQueue(p, ready_queue)&&isReadyQueueSortedByArrivalTime(ready_queue)) {
                ready_queue.push_back(&p);
                
            }
        }
    

        

        // If no process is ready, increment the current time
        if (ready_queue.empty()) {
            ++current_time;
            continue;
        }
        // Get the next process from the ready queue
        Process* p = ready_queue.front();
        ready_queue.erase(ready_queue.begin());

        if (p->start_time == -1) {
            p->start_time = current_time; // Set start time if not already set
        }

        // Determine execution time for the current time slice
        int exec_time = min(p->remaining_time, time_quantum);

        // Run the process for the time slice
        for (int t = current_time; t < current_time + exec_time && t < last_instant; ++t) {
            p->isRunning[t] = 1; // Process is running
        }

        current_time += exec_time;
        p->remaining_time -= exec_time;

        // If the process finishes
        if (p->remaining_time == 0) {
            p->finish_time = current_time;
            p->turnaround = p->finish_time - p->arrival_time;
            p->isCompleted = true;
            ++completed_processes;
        } else {
            // Add the process back to the ready queue if it's not finished
            for (Process& x : processes) {
                if (x.arrival_time <= current_time && !x.isCompleted && !inQueue(x, ready_queue)&& x.name != p->name) {
                    ready_queue.push_back(&x);
                    
                }
            }
            ready_queue.push_back(p);
            
        }
        
    
    }
    // Output results
    if (trace_visualization) {
        trace(processes, last_instant, "RR-" + to_string(time_quantum) );
    } else {
        stats(processes, "RR-"+ to_string(time_quantum));
    }
}

void HRRN(vector<Process> processes, int last_instant, bool trace_visualization) {
    // Sort processes by arrival time initially
    sort(processes.begin(), processes.end(), ArrivalTimeComparator());

    int current_time = 0;
    int completed_processes = 0;

    while (completed_processes < n_of_processes) {
        vector<Process*> ready_queue;

        // Add processes that have arrived by the current time and are not completed
        for (auto& p : processes) {
            if (p.arrival_time <= current_time && !p.isCompleted) {
                ready_queue.push_back(&p);
            }
        }

        // If there are no ready processes, move the time forward
        if (ready_queue.empty()) {
            current_time++;
            continue;
        }

        // Calculate response ratio for each process in the ready queue
        Process* selected_process = nullptr;
        double max_response_ratio = -1.0;

        for (Process* p : ready_queue) {
            double response_ratio = 
                (static_cast<double>(current_time - p->arrival_time) + p->service_time) / p->service_time;

            if (response_ratio > max_response_ratio) {
                max_response_ratio = response_ratio;
                selected_process = p;
            }
        }

        // Execute the selected process
        if (selected_process) {
            selected_process->start_time = (selected_process->start_time == -1) ? current_time : selected_process->start_time;
            current_time += selected_process->service_time;
            selected_process->finish_time = current_time;
            selected_process->turnaround = selected_process->finish_time - selected_process->arrival_time;

            // Update isRunning state
            for (int t = selected_process->arrival_time; t < selected_process->finish_time; ++t) {
                if (t >= selected_process->start_time) {
                    selected_process->isRunning[t] = 1; // Running
                } else {
                    selected_process->isRunning[t] = 0; // Waiting
                }
            }

            // Mark as completed
            selected_process->isCompleted = true;
            completed_processes++;
        }
    }

    // Output results
    if (trace_visualization) {
        trace(processes, last_instant, "HRRN");
    } else {
        stats(processes, "HRRN");
    }
}

void Aging(vector<Process> processes, int last_instant, int quantum, bool trace_visualization) {
    vector<Process*> ready_queue;
    Process* selected_process;
    int current_time = 0;

    // Initialize priorities for processes
    for (auto& process : processes) {
        process.priority = process.service_time ; 
        process.finish_time = last_instant;
    }
    for (auto& process : processes) {
        if (process.arrival_time <= current_time && !process.isCompleted && !inQueue(process, ready_queue)) {
            ready_queue.push_back(&process);
            //process.priority +=1;
        }
    }
    
    while (current_time < last_instant) {
    
        // Admit new processes into the ready queue
        for(int j=0;j<n_of_processes;j++){
            Process *p = &processes[j];
            if (p->arrival_time <= current_time && !inQueue(*p,ready_queue)){
                p->priority++;
                ready_queue.push_back(p);
            }
        }

        if (!ready_queue.empty()) {
            sort(ready_queue.begin(), ready_queue.end(), priorityComparator());

            // Select the process with the highest priority and 
            selected_process = ready_queue.front();
            ready_queue.erase(ready_queue.begin());
            
            // Execute the process for up to 'quantum' time units
            int executed_time = 0;
            while (executed_time < quantum ) {
                selected_process->isRunning[current_time] = 1;// Mark process as running
                //mark the rest of processes as waiting
                updateToWaitQueue(ready_queue, *selected_process, current_time );
                for (auto& process : ready_queue) {
                    if(process != selected_process )
                        process->priority++;
                }
                current_time++;
                executed_time++;
                }
               
                // Reset its priority to the initial value
                selected_process->priority = selected_process->service_time;
                ready_queue.push_back(selected_process);
        } else {
            // If the ready queue is empty, move to the next time step
            current_time++;
        }
    }
    // Visualization or stats output
    if (trace_visualization)
        trace(processes, last_instant, "Aging");
    
}



int main(){
    vector<Process> processes;
    vector <string> schedules;
    // Reading first three inputs
    cin >> visualization;

    cin >> line;
    //schedule = stoi(line);

    // Use a stringstream to parse each line
    stringstream ss(line); //initializes the stringstream object ss with the contents of the string line
    string temp;

    // Extract elements separated by commas
    while (getline(ss, temp, ',')) {
       // cout<<temp<<endl;
        schedules.push_back(temp); // Add parsed value to the vector
    }

    cin >> line;
    last_instant = stoi(line);

    cin >> line;
    n_of_processes = stoi(line);


    // Read lines until EOF
    while (getline(cin, line)) {
        char process_name;
        int arrival_time, service_time;

        if (line.find_first_not_of(" \t\n\r") == string::npos) {
            continue; // Skip this iteration if the line is empty or just whitespace
        }

        // Use a stringstream to parse each line
        stringstream ss(line); //initializes the stringstream object ss with the contents of the string line
        string temp;

        getline(ss, temp, ','); // Extract character
        process_name = temp[0];

        getline(ss, temp, ','); // Extract first integer
        arrival_time = stoi(temp);

        getline(ss, temp, ','); // Extract second integer
        service_time = stoi(temp);

        // Store the parsed data into processed vector
        Process p =Process(process_name, arrival_time, service_time,last_instant);
        processes.push_back(p);
        
    }

    size_t j=0;
    int algorithm,quantum;
    while (j<schedules.size()){


        string schedule = schedules[j];
        //cout<<schedule<<endl;
        if(schedule.size() == 1){
            algorithm = stoi(schedule);
        }
        else{
            stringstream ss(schedule);
            string temp;

            getline(ss, temp, '-');
            algorithm = stoi(temp);

            getline(ss, temp, '-');
            quantum = stoi(temp);  
        }



        switch(algorithm){
            case 1: 
            {
                FCFS(processes, last_instant, visualization == "trace");
                break;
            }
            case 2: 
            { 
                RoundRobin(processes, last_instant, quantum, visualization == "trace");
                break;
            }
            case 3:
            {
                SPN(processes,last_instant, visualization == "trace");
                break;
            }
            case 4:
            {
                SRT(processes,last_instant, visualization == "trace");
                break;
            }
            case 5: 
            {
                HRRN(processes, last_instant, visualization == "trace");
                break;
            }
            case 6:
            {
                FB1(processes,last_instant, visualization == "trace");
                break;
            }
            case 7:
            {
                FB2i(processes,last_instant, visualization == "trace");
                break;
            }
            case 8:
            {
                Aging(processes,last_instant, quantum, visualization == "trace");
                break;
            }

        }
        j++;
    }
        



    return 0;
}
