package TA;

public class LinkedListTA {
	
	private static class Node{
		int data;
		Node next;
		Node(int data){
			this.data = data;
		}
	}
	
	 Node head;
	 Node tail;
	
	 int size;
	
	public void addLast(int data) {
		if(size == 0) {
			Node nn = new Node(data);
			head = tail = nn;
			size++;
		}else {
			Node nn = new Node(data);
			tail.next = nn;
			nn.next = null;
			tail = nn;
			size++;
		}
	}
	
	public void display() {
		Node temp = head;
		while(temp!=null) {
			System.out.print(temp.data+"->");
			temp = temp.next;
		}
	}
	
	public void removeFirst() {
		if(size==0) {
			System.out.println("List is empty");
		}else {
			Node rm = head;
			Node nNode = rm.next;
			rm.next = null;
			head = nNode;
			size--;
		}
	}
	
	public void removeLast() {
		if(size==0) {
			System.out.println("Cannot remove anything");
		}else {
			Node temp = head;
			while(temp.next!=tail) {
				temp = temp.next;
			}
			
			temp.next=null;
			tail = temp;
			size--;
		}
	}
	public void removeAt(int idx) {
		if(idx < 0 || idx>=size) {
			System.out.println("Invalid Index");
		}else if(idx == 0) {
			removeFirst();
		}else if(idx == size-1) {
			removeLast();
		}else {
			Node temp = head;
			for(int i=0; i<idx-1; i++) {
				temp = temp.next;
			}
			
			Node nextNode = temp.next.next;
			temp.next = nextNode;
			size--;
		}
	}
	
	public Node getAt(int idx) {
		Node temp = head;
		
		for(int i=0; i<idx; i++) {
			temp = temp.next;
		}
		
		return temp;
	}
	
	public void reverseData() {
		int left = 0;
		int right = size-1;
		
		while(left<right) {
			Node lnode = getAt(left);
			Node rnode = getAt(right);
			
			int temp = lnode.data;
			lnode.data = rnode.data;
			rnode.data = temp;
			
			left++;
			right--;
		}
	}
	
	public void reversePointer() {
		Node slow=head;
		Node fast=head.next;
		
		while(fast!=null) {
			Node current = fast.next;
			fast.next = slow;
			slow = fast;
			fast = current;
		}
		
		// last me head tail swap karna yad rakhe  , this has to be done explictly, since head and tail are sepeate nodes and static
		Node temp = head;
		head = tail;
		tail = temp;
		
		tail.next = null;
	}
	
	// kth node from last in a single traversal
	public Node kthNode(int k) {
		Node slow = head;
		Node fast = head;
		
		while(k!=0) {
			fast = fast.next;
			k--;
		}
		
		while(fast!=null) {
			slow = slow.next;
			fast = fast.next;
		}
		
		return slow;
		
	}
	
	public int getMid() {
		Node slow = head;
		Node fast = head;
		
		while(fast!=null && fast.next!=null) {
			slow = slow.next;
			fast = fast.next.next;
		}
		
		return slow.data;
	}
	
	public LinkedListTA mergeToSortedLL(Node l1, Node l2) {
		LinkedListTA ll = new LinkedListTA();
		
		Node h1 = l1;
		Node h2 = l2;
		
		while(h1!=null && h2!=null) {
			if(h1.data < h2.data) {
				ll.addLast(h1.data);
				h1 = h1.next;
			}else if(h1.data > h2.data) {
				ll.addLast(h2.data);
				h2 = h2.next;
			}else {
				ll.addLast(h1.data);
				h1= h1.next;
				h2 = h2.next;
			}
		}
	
		while(h1!=null) {
			ll.addLast(h1.data);
			h1=h1.next;
		}
		
		while(h2!=null) {
			ll.addLast(h2.data);
			h2=h2.next;
		}
		
		return ll;
		
	}
	
	public int getFirst() {
		if(this.size == 0) {
			return -1;
		}
		
		return this.head.data;
	}
	
	
	// space complexity remains constant here, because we either add/remove
	public void removeDuplicates() {
		Node nh = null;
		Node nt = null;
		int ns = 0;
		
		while(this.size > 0) {
			int data = this.getFirst();
			
			if(nh == null) {
				Node nn = new Node(data);
				nh = nn;
				nt = nn;
				ns++;
			}else if(nt.data !=data) {
				Node nn = new Node(data);
				nt.next = nn;
				nt = nn;
				ns++;
			}
			
			this.removeFirst();
		}
		
		this.head = nh;
		this.tail = nt;
		this.size = ns;
		
	}
	public static void main(String[] args) {
		
		LinkedListTA l1 = new LinkedListTA();
	
		l1.addLast(10);
		l1.addLast(20);
		l1.addLast(30);
		l1.addLast(40);
		l1.addLast(50);
		l1.addLast(60);
		//addLast(70);
		
		System.out.println("Hello LinkedList");
	
		//reversePointer();
		//display();

//		Node nn = kthNode(5);
//		System.out.println(nn.data);
		
		
		
		LinkedListTA l3 = new LinkedListTA();
		l3.addLast(10);
		l3.addLast(20);
		l3.addLast(30);
		l3.addLast(90);
		
		
		LinkedListTA l2 = new LinkedListTA();
		l2.addLast(2);
		l2.addLast(2);
		l2.addLast(4);
		l2.addLast(4);
		l2.addLast(4);
		l2.addLast(5);
		
	//	System.out.println(l3.head.data + "," +l2.head.data);
		
	//	LinkedListTA ans = l3.mergeToSortedLL(l3.head, l2.head);
		
	//	ans.display();
		
		l2.removeDuplicates();
		l2.display();

	}
	
	

}
