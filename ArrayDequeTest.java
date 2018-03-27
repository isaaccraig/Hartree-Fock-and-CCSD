
public class ArrayDequeTest {

  private ArrayDeque ardeq;

  /** Tests **/
  public ArrayDequeTest() {
      ardeq = new ArrayDeque();
  }

  public boolean TestSize() {
      boolean passed = (ardeq.size() == 0);
      passed = ardeq.isEmpty() && passed;
      ardeq.addFirst(3);
      ardeq.addLast(4);
      passed = (ardeq.size() == 2) && passed;
      passed = (ardeq.get(0) == 3) && (ardeq.get(1) == 4) && passed;
      int first = ardeq.removeFirst();
      passed = (first == 3) && (ardeq.size() == 1) && passed;
      return passed;
  }

  public void main() {

      boolean passed = TestSize();
      if (passed) {
          System.out.println("TEST 1 PASSED !");
      }
      
  }


}
