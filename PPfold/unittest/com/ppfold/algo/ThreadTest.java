package com.ppfold.algo;

import java.util.concurrent.Executor;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicBoolean;

import junit.framework.TestCase;

public class ThreadTest extends TestCase {
    private static final int SIMULTANEOUS = 8;
    
    public void testNewThread() {
        System.out.println(getName());
        Executor executor = new Executor() {
            public void execute(final Runnable command) {
                Thread t = new Thread() {
                    public void run() {
                        command.run();
                    }
                };
                t.start();
            }
        };
        calculate(executor, SIMULTANEOUS);
    }
    
    public void testThreadPool_1() {
        System.out.println(getName());
        Executor executor = new Executor() {
            private ExecutorService threadPool = Executors.newFixedThreadPool(1);
            public void execute(final Runnable command) {
                threadPool.execute(command);
            }
        };
        calculate(executor, SIMULTANEOUS);
    }
    
    public void testThreadPool_2() {
        System.out.println(getName());
        Executor executor = new Executor() {
            private ExecutorService threadPool = Executors.newFixedThreadPool(2);
            public void execute(final Runnable command) {
                threadPool.execute(command);
            }
        };
        calculate(executor, SIMULTANEOUS);
    }
    
    public void testThreadPool_4() {
        System.out.println(getName());
        Executor executor = new Executor() {
            private ExecutorService threadPool = Executors.newFixedThreadPool(4);
            public void execute(final Runnable command) {
                threadPool.execute(command);
            }
        };
        calculate(executor, SIMULTANEOUS);
    }
    
    public void testThreadPool_8() {
        System.out.println(getName());
        Executor executor = new Executor() {
            private ExecutorService threadPool = Executors.newFixedThreadPool(8);
            public void execute(final Runnable command) {
                threadPool.execute(command);
            }
        };
        calculate(executor, SIMULTANEOUS);
    }
    
    public void testSerial() {
        System.out.println(getName());
        Executor executor = new Executor() {
            public void execute(final Runnable command) {
                command.run();
            }
        };
        calculate(executor, SIMULTANEOUS);
    }
    
    void calculate(Executor e, int noOfCommands) {
    	MyRunnable[] commands = new MyRunnable[noOfCommands];
    	for (int i = 0; i<noOfCommands;i++) {
    		commands[i] = createCommand();
    	}
        
        long before = System.currentTimeMillis();
        for (MyRunnable com : commands){
        	e.execute(com);
        }
        
        while(!allFinished(commands)) {
            try {
                Thread.sleep(20);
            }
            catch (InterruptedException e1) {
                e1.printStackTrace();
            }
        }
        System.out.println("Finished in "+(System.currentTimeMillis() - before) + " ms");
    }
    
    private boolean allFinished(MyRunnable[] commands) {
    	for (MyRunnable com : commands){
    		if (!com.isFinished()) {
    			return false;
    		}
    	}
    	return true;
    }
    
    MyRunnable createCommand() {
        final long count = 2000000000L;
        MyRunnable command1 = new MyRunnable() {
            AtomicBoolean isFinished = new AtomicBoolean(false);
            public void run() {
            	JobObjectTest jobobject = new JobObjectTest(2); 
//                System.out.println("run started: "+Thread.currentThread()+ ": "+ System.currentTimeMillis());
                long total = 0;
                for(int c = 0; c<1000000; c++){
                			jobobject.multiplyMatrices();
                		}
//                for (int  i = 0; i< count; i++) {
//                    total += i;
//                }
                isFinished.set(true);
//                System.out.println("run finished: "+Thread.currentThread()+ ": "+ System.currentTimeMillis());
            }
            @Override
            public boolean isFinished() {
                return isFinished.get();
            }
        };
        return command1;
    }

    private static abstract class MyRunnable implements Runnable {
        public abstract boolean isFinished();
    }
}
