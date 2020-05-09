package TA;
import java.util.*;

public class BTreeTA {
	
	// 1. Create Class Node  
	 static class Node{
		int data;
		Node left;
		Node right;
		Node(int data){
			this.data = data;
			this.left = this.right = null;
		}
	}
	
	private static Node root;
	
	// 2. Construct the binary tree
	public static void construct(int[] data) {
		
		ArrayList<Node> cplist = new ArrayList<>();
		
		for(int i=0; i<data.length; i++) {
			if(data[i] == -1) {
				cplist.remove(cplist.size()-1);  
			}else {
				Node nn = new Node(data[i]);
				if(root == null) {
					root = nn;
				}else {
					Node cp = cplist.get(cplist.size()-1);
					if(cp.left == null) {
						cp.left = nn;
					}else {
						cp.right = nn;
					}	
				}
				cplist.add(nn);
			}
		}
	}
	
	// 3. display the binary tree
	public static void display(Node node) {
		
		if(node == null) {
			return;
		}
		
		String s = "";
		
		if(node.left != null) {
			s+=node.left.data + "->";
		}
		
		s+=node.data + "<-";
		
		if(node.right!=null) {
			s+=node.right.data;
		}
		
		System.out.println(s);
		
		if(node.left != null) {
			display(node.left);
		}
		
		if(node.right !=null) {
			display(node.right);
		}
	
	}
	
	// 4. Size
	public static int size(Node node) {
		if(node == null) {
			return 0;
		}
		
		int s = 0;
		
		if(node.left!= null) {
			s+=size(node.left);
		}
		
		if(node.right!=null) {
			s+=size(node.right);
		}
		
		return s+1;
	}
	
	// 5. height
	public static int height(Node node) {
	
		if(node == null) {
			return 0;
		}
		
		int hl = height(node.left);
		
		int rl = height(node.right);
		
		return Math.max(hl, rl) + 1;
		
	}
	
	// 6. Maximum of Tree
	public static int max(Node node) {
		if(node == null) {
			return 0;
		}
		
		int m1 = max(node.left);
		
		int m2 = max(node.right);
		
		return Math.max(Math.max(m1, m2), node.data);
		
	}
	
	// 7. Find the given element in tree and give true/false
	public static boolean find(Node node, int dtf) {
		if(node.data == dtf) {
			return true;
		}
		
		if(node.left!=null) {
			boolean bl = find(node.left,dtf);
			if(bl) {
				return true;
			}
		}
		
		if(node.right!=null) {
			boolean br = find(node.right, dtf);
			if(br) {
				return true;
			}
		}
		
		return false;
	}
	
	// 8. Find the path from the given node to root in the form of ArrayList<Integer>
	public static ArrayList<Integer> n2rpath(Node node, int dtf){
		
		if(node == null) {
			ArrayList<Integer> br = new ArrayList<>();
			return br;
		}
		
		if(node.data == dtf) {
			ArrayList<Integer> nr = new ArrayList<>();
			nr.add(node.data);
			return nr;
		}
		
		if(node.left!=null) {
			ArrayList<Integer> clpath = n2rpath(node.left, dtf);
			if(clpath!=null) {
				clpath.add(node.data);
				return clpath;
			}
		}
		
		if(node.right!=null) {
			ArrayList<Integer> crpath = n2rpath(node.right, dtf);
			if(crpath!=null) {
				crpath.add(node.data);
				return crpath;
			}
		}
		
		return null;
		
	}

	// 9. Traversals : Pre, Post, In
	
	public static void preOrder(Node node) {
		if(node == null) {
			return;
		}
		System.out.print(node.data+" ");
		
		if(node.left!=null) {
			preOrder(node.left);
		}
		
		if(node.right!=null) {
			preOrder(node.right);
		}
	}
	
	public static void postOrder(Node node) {
		if(node == null) {
			return;
		}
		
		if(node.left!=null) {
			postOrder(node.left);
		}
		
		if(node.right!=null) {
			postOrder(node.right);
		}
		
		System.out.print(node.data+" ");
	}
	
	public static void inorder(Node node) {
		if(node == null) {
			return;
		}
		
		if(node.left!=null) {
			inorder(node.left);
		}
		
		System.out.print(node.data+" ");
		
		if(node.right!=null) {
			inorder(node.right);
		}
		
	}
	
	// 10. Construct tree from inorder and preorder
	
	public static Node constructPreIn(ArrayList<Integer> pre, ArrayList<Integer> in, int plow, int phigh, int inlow, int ihigh) {
		
		if(plow>phigh || inlow>ihigh) {
			return null;
		}
		
		Node root = new Node(pre.get(plow));
		
		int lhs = 0;
		while(in.get(lhs+inlow)!=pre.get(plow)) {
			lhs++;
		}
		
		
		root.left = constructPreIn(pre,in,plow+1,plow+lhs,inlow,inlow+lhs-1);
		root.right = constructPreIn(pre,in,plow+lhs+1,phigh,inlow+lhs+1,ihigh);
		
		return root;
	}
	
	
	// pre-order iterative
	
	// 1. Add root to stack
	// 2. Remove the last element
	// 3. Print
	// 4. Add children in reverse order
	
	public static void preOrderIterative(Node node) {
		ArrayList<Node> stack = new ArrayList<>();
		stack.add(node);
		
		while(stack.size() > 0) {
			Node rm = stack.remove(stack.size()-1);
			
			System.out.print(rm.data+" ");
			
			if(rm.right!=null) {
				stack.add(rm.right);
			}
			
			if(rm.left!=null) {
				stack.add(rm.left);
			}
			
			
		}
	}
	
	// 1. O(n*n) complexity
	
	public static int diameter(Node node) {
		
		if(node == null) {
			return 0;
		}
		
		int lHeight = height(node.left);
		int rHeight = height(node.right);
		
		int lMax = diameter(node.left);
		int rMax = diameter(node.right);
		
		return Math.max(lHeight+rHeight + 1, Math.max(lMax, rMax));
	}
	
	// 2. O(n) complexity
	
	private static class DiaPair{
		Node n=null;
		int ht=0;
		int dia=0;
		boolean f=true;
	}

	public static DiaPair diameter2(Node root) {
		
		if(root==null) {
			DiaPair bp=new DiaPair();
			bp.ht=0;
			bp.dia=0;
			return bp;
		}
		
		DiaPair l=diameter2(root.left);
		DiaPair r=diameter2(root.right);
		
		DiaPair mr=new DiaPair();
		mr.ht=Math.max(l.ht, r.ht)+1;
		mr.dia=Math.max(l.ht+r.ht+1, Math.max(l.dia, r.dia));
		return mr;
	}
	
	public static DiaPair isBalanced(Node node) {
		if(node == null) {
			DiaPair bp = new DiaPair();
			return bp;
		}
		
		DiaPair dl = isBalanced(node.left);
		DiaPair dr = isBalanced(node.right);
		
		DiaPair mp = new DiaPair();
		
		mp.ht = Math.max(dl.ht, dr.ht)+1;
		
		int param = dl.ht - dr.ht;
		
		if(Math.abs(param)>1 && !dl.f && !dr.f) {
			mp.f = false;
		}
		
		return mp;
	}
	
	public static class BstPair{
		int min;
		int max;
		boolean f=true;
		int bstData;
		int bstSize;
	}
	public static BstPair isBST(Node node) {
		if(node == null) {
			BstPair bp = new BstPair();
			bp.min = Integer.MAX_VALUE;
			bp.max = Integer.MIN_VALUE;
			return bp;
		}
		
		BstPair mp = new BstPair();
		
		BstPair lp = isBST(node.left);
		BstPair rp = isBST(node.right);
		
		mp.max = Math.max(lp.max, Math.max(rp.max, node.data));
		
		mp.min = Math.min(lp.min, Math.min(rp.min, node.data));
		
		mp.f = (node.data > lp.max) && (node.data < rp.min) && lp.f && rp.f;
		
		return mp;
	}
	
	
	public static BstPair largestBST(Node node) {
		if(node == null) {
			BstPair bp = new BstPair();
			bp.bstSize = 0;
			bp.bstData = -1;
			bp.max = Integer.MIN_VALUE;
			bp.min = Integer.MAX_VALUE;
			return bp;
		}
		
		BstPair lp = largestBST(node.left);
		BstPair rp = largestBST(node.right);
		
		BstPair mp = new BstPair();
		
		mp.min = Math.min(node.data, Math.min(lp.min, rp.min));
		mp.max = Math.max(node.data, Math.max(lp.max, rp.max));
		
		mp.f = (node.data > lp.max) && (node.data < rp.min) && lp.f && rp.f;
		
		if(mp.f == true) {
			mp.bstSize = lp.bstSize + rp.bstSize +1;
			mp.bstData = node.data;
		}else {
			if(lp.bstSize > rp.bstSize) {
				mp.bstSize = lp.bstSize;
				mp.bstData = lp.bstData;
			}else {
				mp.bstSize = rp.bstSize;
				mp.bstData = rp.bstData;
			}
		}
		
		return mp;
	}

	public static int countLeaves(Node node) {
		if(node == null) {
			return 0;
		}
		
		int cnt = 0;
		
		if(node.left==null && node.right == null) {
			cnt++;
		}
		
		cnt+=countLeaves(node.left);
		cnt+=countLeaves(node.right);
		
		return cnt;
	}
	
	 public static int sumOfLeftLeaves(Node root) {
	       if(root == null) {
	    	   return 0;
	       }
	       
	       int sum = 0;
	       if(root.left!=null) {
	    	   sum+=root.left.data;
	       }
	       
	       
	       sum+=sumOfLeftLeaves(root.left);
	       sum+=sumOfLeftLeaves(root.right);
	       
	       return sum;
	       
	 }
	 
	 public static boolean isLeafAtSameLevel(Node root) {
		 
		 if(root == null) {
			 return true;
		}
		 
		return checkHelper(root,level);
	 }
	
	 static int level = -1;
	 
	private static boolean checkHelper(Node root, int l) {
		
		if(root == null) {
			return true;
		}
		
		// check if leaf nodes are at the same level
		if(root.left == null && root.right == null) {
			if(level == -1) {
				level = l;
				return true;
			}else {
				if(level != l) {
					return false;
				}
				
				return true;
			}
		}
		
		// bring the result from left node
	
		boolean lres = checkHelper(root.left, l+1);
		
		if(lres == false) {
			return false;
		}
		
		// bring the result from right node
		boolean rres = checkHelper(root.right, l+1);
		
		if(rres == false) {
			return false;
		}
		
		return true;
		
	}
	
	public static int countNonLeafNodes(Node node) {
		if(node == null) {
			return 0;
		}
		
		if(node.left == null && node.right == null) {
			return 0;
		}
		
		
		int lres = countNonLeafNodes(node.left);
		
		int rres = countNonLeafNodes(node.right);
		
		
		return lres + rres + 1;
	}
	
	// get all the non leaf nodes of a tree in a seperate arrayList
	
	public static ArrayList<Integer> nonLeafNodes(Node node){
		if(node == null) {
			return new ArrayList<Integer>();
		}
		
		ArrayList<Integer> myres = new ArrayList<Integer>();
		
		if(node.left!=null || node.right!=null) {
			myres.add(node.data);
		}
		
		ArrayList<Integer> lRes = nonLeafNodes(node.left);
		
		ArrayList<Integer> rRes = nonLeafNodes(node.right);
		
		myres.addAll(lRes);
		myres.addAll(rRes);
		
		return myres;	
	}
	
	// get all root to leaf paths in binary tree
	
	public static ArrayList<String> node2LeafPaths(Node root){
		ArrayList<String> list = new ArrayList<String>();
		return node2LeafHelper(root,"" ,list);
	}
	
	// we are solving this question using pre-order traversal
	private static ArrayList<String> node2LeafHelper (Node node, String psf, ArrayList<String> ls){
		if(node == null) {
			return ls;
		}
		
		// leaf node is the closing case here , i.e , base case 1
		// in the case of leaf node, just add that node to path and return arraylist
		
		if(node.left == null && node.right==null) {
			psf += node.data;
			ls.add(psf);
			return ls;
		}
		
		psf  = psf + node.data +"->";
		
		ls = node2LeafHelper(node.left, psf, ls);
		ls = node2LeafHelper(node.right, psf, ls);
		
		return ls;
	}
	
	
	// LCA O(n)
	
	static int lca = 0;
	
	public static Node LCA(Node root, int n1, int n2) {
		
		if(root == null) {
			return null;
		}
		
		if(root.data == n1 || root.data == n2) {
			return root;
		}
		
		Node lans = LCA(root.left, n1, n2);
		Node rans = LCA(root.right, n1, n2);
		
		if(lans !=null && rans!=null) {
			return root;
		}
		
		if(lans!=null) {
			return lans;
		}else {
			return rans;
		}
	}
	
	public static void levelOrder(Node root) {
		Queue<Node> qu = new LinkedList<>();
		qu.add(root);
		
		ArrayList<Integer> lo = new ArrayList<Integer>();
		
		while(qu.size()>0) {
			
			Node rm = qu.remove();
			
			if(rm!=null) {
				lo.add(rm.data);
			}
			
			if(rm.left!=null) {
				qu.add(rm.left);
			}
			
			if(rm.right!=null) {
				qu.add(rm.right);
			}
		}
		
		System.out.println(lo);
	}
	
	// **************************************************************************************************************//
	// Binary Tree Views    -> Always try level order traversal
	
	// 1. Right Side View
	
	
	public static ArrayList<Integer> rightView(Node node) {
		Queue<Node> qu = new LinkedList<Node>();
		qu.add(node);
		qu.add(null);
		
		ArrayList<Integer> lo = new ArrayList<Integer>();
		ArrayList<Integer> rw = new ArrayList<Integer>();
		
		while(qu.size()>0) {
			Node curr = qu.remove();
			
			if(curr!=null) {
				lo.add(curr.data);
				if(curr.left!=null) {
					qu.add(curr.left);
				}
				
				if(curr.right!=null) {
					qu.add(curr.right);
				}
				
			}else {
				if(qu.size()>0 && lo.size()>0) {	
					int num = lo.remove(lo.size()-1);
					rw.add(num);
					qu.add(null);
				}
			}
		}
		
		if(lo.size()>0) {
			rw.add(lo.get(lo.size()-1));
		}
		
		return rw;
	}
	
	// 2. Left View
	public static ArrayList<Integer> leftView(Node node) {
		Queue<Node> qu = new LinkedList<Node>();
		qu.add(node);
		qu.add(null);
		
		ArrayList<Integer> lo = new ArrayList<Integer>();
		ArrayList<Integer> lw = new ArrayList<Integer>();
		
		lw.add(node.data);
		
		while(qu.size()>0) {
			Node curr = qu.remove();
			
			if(curr!=null) {
				lo.add(curr.data);
				if(curr.left!=null) {
					qu.add(curr.left);
				}
				
				if(curr.right!=null) {
					qu.add(curr.right);
				}
				
			}else {
				if(qu.size()>0){
					lo.add(-1);
					qu.add(null);
				}
			}
		}
		
		for(int i=0; i<lo.size(); i++) {
			int ele = lo.get(i);
			if(ele == -1 && i<lo.size()) {
				lw.add(lo.get(i+1));
			}
		}
		
		return lw;
		
	}
	
	// 3. Top View , using vertical order
	// add karte rho, already present h hashmap me to kuchmat karo bas children add karo
	
	public static class TPair{
		Node node;
		int level;
		
		TPair(Node node, int level){
			this.node = node;
			this.level = level;
		}
	}
	
	public static void topView(Node root) {
		Queue<TPair> qu = new LinkedList<>();
		HashMap<Integer,Integer> hmap = new HashMap<>();
		
		TPair rp = new TPair(root,0);
		qu.add(rp);
		
		while(qu.size()>0) {
			TPair rm = qu.remove();                  // remove
			
			if(!hmap.containsKey(rm.level)) {        // check hashmap 
				hmap.put(rm.level, rm.node.data);
			}
			
			if(rm.node.left!=null) {
				TPair np = new TPair(rm.node.left, rm.level-1);
				qu.add(np);
			}
			
			if(rm.node.right!=null) {
				TPair np = new TPair(rm.node.right, rm.level+1);
				qu.add(np);
			}
			
		}
		
		
		// now sort this hmap on the basis of keyvalue
		
		TreeMap<Integer, Integer> sorted = new TreeMap<>(); 
		 sorted.putAll(hmap); 
		 
		 for (Map.Entry<Integer, Integer> entry : sorted.entrySet())  
	            System.out.print(entry.getValue()+" ");  
		 
		 
	}
	
	// add karte rho, already present h hashmap me to update karte rho aur children add karo
	
	public static void bottomView(Node root) {
		Queue<TPair> qu = new LinkedList<>();
		HashMap<Integer,Integer> hmap = new HashMap<>();
		
		TPair rp = new TPair(root,0);
		qu.add(rp);
		
		while(qu.size()>0) {
			TPair rm = qu.remove();                  // remove
			
			if(!hmap.containsKey(rm.level)) {        // check hashmap 
				hmap.put(rm.level, rm.node.data);
			}else {
				hmap.put(rm.level, rm.node.data);
			}
			
			if(rm.node.left!=null) {
				TPair np = new TPair(rm.node.left, rm.level-1);
				qu.add(np);
			}
			
			if(rm.node.right!=null) {
				TPair np = new TPair(rm.node.right, rm.level+1);
				qu.add(np);
			}
			
		}
		
		
		// now sort this hmap on the basis of keyvalue
		
		TreeMap<Integer, Integer> sorted = new TreeMap<>(); 
		 sorted.putAll(hmap); 
		 
		 for (Map.Entry<Integer, Integer> entry : sorted.entrySet())  
	            System.out.print(entry.getValue()+" ");  
	}
	
	// ik hashmap arraylist and levels ka banao, sab add karte rho.
	
	public static void verticalView(Node root) {
		Queue<TPair> qu = new LinkedList<>();
		HashMap<Integer,ArrayList<Integer>> hmap = new HashMap<>();
		
		TPair rp = new TPair(root,0);
		qu.add(rp);
		
		while(qu.size()>0) {
			TPair rm = qu.remove();                  // remove
			
			
			if(hmap.containsKey(rm.level)) {
				ArrayList<Integer> list = hmap.get(rm.level);
				list.add(rm.node.data);
			}else {
				ArrayList<Integer> nlist = new ArrayList<>();
				nlist.add(rm.node.data);
				hmap.put(rm.level, nlist);
			}
			
			
			if(rm.node.left!=null) {
				TPair np = new TPair(rm.node.left, rm.level-1);
				qu.add(np);
			}
			
			if(rm.node.right!=null) {
				TPair np = new TPair(rm.node.right, rm.level+1);
				qu.add(np);
			}
			
		}
		
		
		// now sort this hmap on the basis of keyvalue
		
		ArrayList<ArrayList<Integer>> ans = new ArrayList<ArrayList<Integer>>();
		
		TreeMap<Integer, ArrayList<Integer>> sorted = new TreeMap<>(); 
		 sorted.putAll(hmap); 
		 
		 for (Map.Entry<Integer, ArrayList<Integer>> entry : sorted.entrySet())  {
			 ArrayList<Integer> nlist = entry.getValue();
			 ans.add(nlist);
		 }
		 
		 
		System.out.println(ans);

	}
	
	// left child-> diag +1, right child-> with same diagonal value
	
	// this is according to geeksforgeeks
	
	public static class DiagPair{
		Node node;
		int diag;
		
		DiagPair(Node node, int diag){
			this.node = node;
			this.diag = diag;
		}
	}
	
	public static void diagonalTraversal(Node root) {
		Queue<DiagPair> qu = new LinkedList<>();
		
		HashMap<Integer, ArrayList<Integer>> hmap = new HashMap<>();
		
		DiagPair rp = new DiagPair(root, 0);
		qu.add(rp);
		
		while(qu.size() > 0) {
			DiagPair rm = qu.remove();
			
			if(hmap.containsKey(rm.diag)) {
				ArrayList<Integer> clist = hmap.get(rm.diag);
				clist.add(rm.node.data);
				hmap.put(rm.diag, clist);
			}else {
				ArrayList<Integer> nlist = new ArrayList<Integer>();
				nlist.add(rm.node.data);
				hmap.put(rm.diag, nlist);
			}
			
			if(rm.node.left!=null) {
				DiagPair np = new DiagPair(rm.node.left, rm.diag-1);
				qu.add(np);
			}
			
			if(rm.node.right!=null) {
				DiagPair np = new DiagPair(rm.node.right, rm.diag);
				qu.add(np);
			}
			
		}
		
		
        ArrayList<ArrayList<Integer>> ans = new ArrayList<ArrayList<Integer>>();
		
		TreeMap<Integer, ArrayList<Integer>> sorted = new TreeMap<>(); 
		 sorted.putAll(hmap); 
		 
		 for (Map.Entry<Integer, ArrayList<Integer>> entry : sorted.entrySet())  {
			 ArrayList<Integer> nlist = entry.getValue();
			 ans.add(nlist);
		 }

		
		for(int i=0; i<ans.size(); i++) {
			ArrayList<Integer> list = ans.get(i);
	
			for(int ele: list) {
				System.out.print(ele+" ");
			}
			System.out.println();
		}
		
		//System.out.println(ans);
	}
	
	static ArrayList<Integer> ans = new ArrayList<Integer>();
	
	public static void boundaryTraversal(Node root) {
		
		if(root == null) {
			return ;
		}
		
		ans.add(root.data);
		
		// 1. call for printing the left boundary
		
		if(root.left!=null) {
			leftBoundary(root.left);
		}
		
		// 2. call for printing the leaf nodes
		
		printLeafs(root.left);
		
		printLeafs(root.right);
		
		// 3. call for printing the right boundary, in reverse order
		
		if(root.right!=null) {
			rightBoundary(root.right);
		}
	
	}
	
	
	public static void leftBoundary(Node node) {
		if((node == null) || (node.left == null && node.right==null)) {
			return ;
		}
	
		ans.add(node.data);
		
		if(node.left==null) {
			leftBoundary(node.right);
		}else {
			leftBoundary(node.left);
		}
		
	}
	
	public static void rightBoundary(Node node) {
		if((node == null) || (node.left== null && node.right == null)) {
			return;
		}
		
		if(node.left == null) {
			rightBoundary(node.right);
		}else {
			rightBoundary(node.left);
		}
		
	   ans.add(node.data);
	
	}
	
	
	public static void printLeafs(Node node) {
		if(node == null) {
			return;
		}
		
		if(node.left==null && node.right==null) {
			ans.add(node.data);
		}
		
		printLeafs(node.left);
		printLeafs(node.right);
		
		return ;
	}
	
	
	// distribute the coins in btree such that each node has 1 coin, and return such moves for distribution
	
	 static int moves = 0;
     public static int distributeCoins(Node root) {
        if(root==null) {
        	return 0;
        }
        
        int lans = 0;
        int rans = 0;
        
        if(root.left!=null) {
           lans = distributeCoins(root.left);
        }
        
        if(root.right!=null) {
        	rans = distributeCoins(root.right);
        }
        
        moves = Math.abs(lans) + Math.abs(rans);
        
        return lans+rans+root.data-1;
    }
	
     
     public static String serialize(Node root) {
    	 if(root==null) {
    		 return "";
    	 }
    	 
    	 String ans = "";
    	 
    	ans = helper1(root, ans);
    	 
    	return ans;
     }
     
     public static String helper1(Node root, String ans) {
    	if(root==null) {
    		return "#";
    	}
    	
    	String lans = helper1(root.left,ans);
    	String rans = helper1(root.right,ans);
    	
    	ans = root.data+lans+rans;
    	
    	return ans;
     }
     
     public static Node deserialize(String path) {
    	 if(path.length()==0) {
    		 return null;
    	 }
    	 
    	 Node root ;
    	 
    	 root = helper2(path);
    	 
    	 return root;
     }
   
   
	static int idx = 0;
	private static Node helper2(String path) {
		
		if(idx==path.length()-1) {
			return root;
		}
		
		char dta = path.charAt(idx);
		if(dta=='#') {
			idx=+2;
			return null;
		}else {
			Node lnode = helper2(path.substring(idx+1));
			Node rnode = helper2(path.substring(idx+1));
			
			Node node = new Node(dta-'0');
			node.left = lnode;
			node.right = rnode;
			
			return node;
		}
		
	}
	
	// find duplicate subtrees in a given binary tree
	// using hasmap, strore all the strings while serializing the tree
	

	public static void main(String[] args) {
		//int[] data = {10,20,30,80,-1,-1,40,-1,-1,50,70,-1,60};
		
		//int[] data2 = {1,2,4,8,-1,-1,5,-1,3,6,-1,7};
		
		int[] data1 = {3,2,1,-1,7,-1,-1,5,2,-1,8};
		construct(data1);
		//display(root);
		
		//System.out.println(size(root));
		
		//System.out.println(height(root));
		
		//System.out.println(max(root));
		
		//System.out.println(find(root,60));
		
		//System.out.println(n2rpath(root,60));
		
		//preOrder(root);
		
		//postOrder(root);
		
		//preOrderIterative(root);
		
		//BstPair res = largestBST(root);
		
		//System.out.println(res.bstSize+ " "+res.bstData);

		//System.out.println(countLeaves(root));
		
		//System.out.println(isLeafAtSameLevel(root));
		
		//System.out.println(countNonLeafNodes(root));
		
		//System.out.println(nonLeafNodes(root));
		
		//System.out.println(node2LeafPaths(root));

		//Node ans = LCA(root, 12, 87);
		
		//System.out.println(ans.data);
		
		//levelOrder(root);
		
		//System.out.println(rightView(root));
		
		//System.out.println(leftView(root));
		
		

		//verticalView(root);
		
		//boundaryTraversal(root);
		
		//System.out.println(ans);
		
		//System.out.println();
		
		String str = serialize(root);
		System.out.println(str);
		Node node = deserialize(str);
		
		System.out.println(node.data);
		
		display(node);
	}
}
