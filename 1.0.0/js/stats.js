document.addEventListener('DOMContentLoaded', function() {
  // Update view counts or other analytics
  // This is a placeholder for future stats functionality
  console.log('Documentation viewed');
  
  // Example function to count tools by category
  function countToolsByCategory() {
    const categories = document.querySelectorAll('h2');
    let counts = {};
    
    categories.forEach(cat => {
      const title = cat.textContent;
      const toolsList = cat.nextElementSibling?.nextElementSibling;
      if (toolsList && toolsList.tagName === 'UL') {
        counts[title] = toolsList.querySelectorAll('li').length;
      }
    });
    
    return counts;
  }
  
  // Add event listeners for navigation tracking if needed
});