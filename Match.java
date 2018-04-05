
import java.util.*;
import java.io.*;	


public class Match 
{
	
	public static void main(String[] arg) 
	{
		
		double startTime = System.currentTimeMillis();
		Runtime runtime = Runtime.getRuntime();	
		String corpusSeq = null;
		
		String patternSeq = null;

		Scanner sc = new Scanner(System.in);
		System.out.println("Enter the path of corpus file:");
		String corpus = sc.nextLine();
		System.out.println("Enter the path of pattern file:");
		String pattern = sc.nextLine();								
		corpusSeq = SeqReader.readSeq(corpus);

		patternSeq = SeqReader.readSeq(pattern);
		
		System.out.println("CORPUS: " + corpusSeq.length() + " bases");

		System.out.println("PATTERN: " + patternSeq.length() + " bases");

		
		System.out.print("Match length? ");
		
		int matchLength = sc.nextInt();

		

		SimplifiedOpenAddressingBucketMapping<String,Integer> mapping =
	new SimplifiedOpenAddressingBucketMapping<String,Integer>(patternSeq.length());


		
for (int j = 0; j < patternSeq.length() - matchLength + 1; j++) 
		{
			
			String tag = patternSeq.substring(j, j + matchLength);

			mapping.put(tag, new Integer(j));

		}
		
		System.out.println("\nAfter creating the table, it holds "
 + mapping.getSize() + " sequences of length " + matchLength);

		System.out.println("The number of distinct substrings is " + mapping.getNumTags());

			
		System.out.println("\nThe matches are:");

		int numMatches = 0;

		int numTagMatches = 0;
		for (int j = 0; j < corpusSeq.length() - matchLength + 1; j++)
		{
			
			String tag = corpusSeq.substring(j, j + matchLength);
			
			Iterator loc = mapping.get(tag);

			if (loc != null)
			{
				
				numTagMatches++;
				
				while (loc.hasNext())
				{
					
					System.out.println(j + " " + loc.next() + " " + tag);
					
					numMatches++;
				
				}
			
			}

		}
		
		System.out.print("\nThere were " + numMatches + " matches found ");

		System.out.println("among " + numTagMatches + " matching substrings.");

		
		long kb = 1024;
		double stopTime = System.currentTimeMillis();
        	double elapsedTime = stopTime - startTime;
        	System.out.println("Running time :" + elapsedTime/1000 + " seconds.");

		long usedMemory = (runtime.totalMemory() - runtime.freeMemory()) / kb;
		System.out.println("Used Memory :" + usedMemory + " KB");		
	}	

}

class SeqReader 
{
    
	public static String readSeq(String fileName) 
	{
	
		BufferedReader r;
	
	
	
	try 
		{
            
			InputStream is = new FileInputStream(fileName);
	
	    	r = new BufferedReader(new InputStreamReader(is));
	
		}
	
		catch (IOException e) 
		{
	    
			System.out.println("IOException while opening " +
 fileName + "\n" + e);
	
	    	return null;
	
		}
	
		
		StringWriter buffer = new StringWriter();
	
	
	
	
		try 
		{
	   
			boolean stop = false;
	    
	    
			while (!stop)
		
			{
		    
				String nextline = r.readLine();
		    
				if (nextline == null)
					stop = true;
		    
				else
			
				{
			    
					String seq = nextline.trim();
		
			buffer.write(seq.toLowerCase());
	
			}
		
			}
	
		}
	
		catch (IOException e) 
		{
	    
			System.out.println("IOException while reading sequence from " +
 fileName + "\n" + e);
		
	return null;
	
		}
	
	
		return buffer.toString();
    
	}

}


class TaggedElement<T,E> {
	T tag;
	ArrayList<E> positions;

	TaggedElement(T t){
		tag = t;
		positions = new ArrayList<E>(1);
	}
}

//A limited version of an Open Addressing Bucket Mapping

class SimplifiedOpenAddressingBucketMapping<T,E> {
	
	public static final Object EMPTY = new EmptySlot(); //sentinel for empty slot
	static class EmptySlot{}
	
	public static final Object DELETED = new DeletedElement(); //deleted sentinel
	static class DeletedElement{}

	public static final double DEFAULT_LOAD = 0.5; //default value for target load
	public static final int DEFAULT_CAPACITY = 8; //default capacity
	double targetLoad; //the desired value for size/table.length
	
	Object[] table;   //the underlying array for the hash table
	int size = 0;     //the number of elements put in the collection
	int d = 0;        //the number of slots marked as deleted
	int numTags = 0;  //the number of tags in the collection
	int minCapacity;  // min capacity from parameter given in the constructor
	int highWaterMark;// max value for number of tags before increasing table size

//Constructors 
	
	public SimplifiedOpenAddressingBucketMapping(int capacity, double load){
		if (load >= 1.0)
			throw new IllegalArgumentException("Load must be < 1");
		this.targetLoad = load;
		int tableSize =(int) Math.pow(2,Math.ceil(Math.log(capacity/load)/Math.log(2)));
		table = new Object[tableSize];
		Arrays.fill(table, EMPTY);
		minCapacity = capacity;
		highWaterMark = (int) (table.length*(1+targetLoad)/2);
	}
	
	public SimplifiedOpenAddressingBucketMapping(){
		this(DEFAULT_CAPACITY, DEFAULT_LOAD);
	}

	public SimplifiedOpenAddressingBucketMapping(int capacity) {
		this(capacity, DEFAULT_LOAD);
	}

// Accessors to return the number of tags or number of elements put into the collection
	public int getNumTags() {
		return numTags;
	}
	
// Returns the number of elements that have been inserted into the bucket mapping
	public int getSize() {
		return size;
	}

// Returns true exactly when slot is in use in the hash table
	boolean inUse(int slot){
		return table[slot] != EMPTY && table[slot] != DELETED;
	}

// primary hash function
	static double A = (Math.sqrt(5.0)-1)/2; //multiplier for hash function
	protected int hash(int hashCode) {
		double frac = (hashCode * A) - (int) (hashCode * A);  //fractional part of x*A
		int hashValue = (int) (table.length * frac);  //multiply by m
		if (hashValue < 0) // if this is negative add m to get
			hashValue += table.length;  // hashValue mod m
		return hashValue;
	}

//secondary hash function
	protected int stepHash(int hashCode) {
		int s = (hashCode % (table.length/2 - 1));
		if (s < 0)
			s += (table.length/2 - 1);
		return 2*s + 1;
	}

// internal method to locate the slot where a given tag is held or where it should be
//   inserted if the tag is not currently used by an element in the collection
	
	final int NONE = -1;
	@SuppressWarnings("unchecked")
	protected int locate(T target) {
		int hashCode = target.hashCode();
		int index = hash(hashCode);      //first slot in probe sequence
		int step = stepHash(hashCode);   //step size between probes
		int insertPosition = NONE;
		while (table[index] != EMPTY) {  //continue until an empty slot found
			if (table[index] != DELETED) {
				if (((TaggedElement<T,E>) table[index]).tag.equals(target))
					return index;
			} else if (insertPosition == NONE) {
				insertPosition = index;
			}
			index = (index + step) % table.length;  //move forward in probe sequence
		}
		if (insertPosition == NONE)
			insertPosition = index;
		return insertPosition;
	}

// Returns true exactly when some element has the given tag.
	public boolean containsTag(T tag) {
		return inUse(locate(tag));
	}

// Returns an iterator over the bucket associated with the given tag (or null if
//   there is no element with this tag).
	@SuppressWarnings("unchecked")
	public Iterator<E> get(T tag) {
		int slot = locate(tag);
		if (inUse(slot))
			return ((TaggedElement<T,E>) table[slot]).positions.iterator(); 
		return null;
	}

// Puts a new tagged element into the bucket mapping with the given tag and element.
	@SuppressWarnings("unchecked")
	public void put(T tag, E element) {
		size++;
		int slot = locate(tag);  // get insert position
		if (!inUse(slot)) { // tag needs to be added
			numTags++;
			table[slot] = new TaggedElement<T,E>(tag);
		}
		((TaggedElement<T,E>) table[slot]).positions.add(element);	
		growTableAsNeeded();   
	}

//Methods used to delete elements	
	boolean clearSlot(int slot) {
		if (inUse(slot)) {
			table[slot] = DELETED;  
			d++;                     
			size--;	                  
			return true;
		} else
			return false;
	}
	
	public boolean removeBucket(T tag) {
		return clearSlot(locate(tag)); 
	}
	
// Methods to handle resizing the table to ensure that the true load isn't too high

	void growTableAsNeeded() {
		if (numTags > highWaterMark)
			resizeTable(Math.max(numTags, minCapacity));		
	}

	@SuppressWarnings("unchecked")
	void resizeTable(int desiredCapacity){
		SimplifiedOpenAddressingBucketMapping<T,E> newTable = 
			new SimplifiedOpenAddressingBucketMapping<T,E>(desiredCapacity, targetLoad);
		for (int slot = 0; slot < table.length; slot++)    //insert all elements
			if (inUse(slot)) {
				TaggedElement<T,E> te = (TaggedElement<T,E>) table[slot];
				int newSlot = newTable.locate((T) te.tag);  // get insert position
				newTable.table[newSlot] = te;
			}
		this.table = newTable.table;
		d = 0;
		highWaterMark = (int) (table.length*(1+targetLoad)/2);
	}
}


