/**
 * Wizard.
 *
 * @author Htmlstream
 * @version 1.0
 *
 */
;(function ($) {
  'use strict';

  $.HSCore.components.HSWizard = {
    /**
     *
     *
     * @var Object _baseConfig
     */
    _baseConfig: {},

    /**
     *
     *
     * @var jQuery pageCollection
     */
    pageCollection: $(),

    /**
     * Initialization of Wizard.
     *
     * @param String selector (optional)
     * @param Object config (optional)
     *
     * @return jQuery pageCollection - collection of initialized items.
     */

    init: function (selector, config) {

      this.collection = selector && $(selector).length ? $(selector) : $();
      if (!$(selector).length) return;

      this.config = config && $.isPlainObject(config) ?
        $.extend({}, this._baseConfig, config) : this._baseConfig;

      this.config.itemSelector = selector;

      this.initWizard();

      return this.pageCollection;

    },

    initWizard: function () {
      //Variables
      var $self = this,
        collection = $self.pageCollection;

      //Actions
      this.collection.each(function (i, el) {
        //Variables
        var $this = $(el),
          wizardSteps = $this.data('wizard-steps'),
          $wizardStepsItems = $(wizardSteps).find('> *'),

          wizardContent = $this.data('wizard-content'),
          $wizardContentItems = $(wizardContent).find('> *'),
          $stepsActiveItem = $(wizardContent).find('> .active');

        $wizardContentItems.not('.active').hide();
        $wizardStepsItems.eq($stepsActiveItem.index()).addClass('active');

        $('[data-next-step]').on('click', function (e) {
          e.preventDefault();

          var $this = $(this),
            nextID = $this.data('next-step');

          $wizardStepsItems.removeClass('active');
          $wizardStepsItems.eq($(nextID).index() - 1).addClass('u-checked');
          $wizardStepsItems.eq($(nextID).index()).addClass('active');

          $wizardContentItems.hide().removeClass('active');
          $(nextID).fadeIn(400).addClass('active');
        });

        $('[data-previous-step]').on('click', function (e) {
          e.preventDefault();

          var $this = $(this),
            prevID = $this.data('previous-step');

          $wizardStepsItems.removeClass('active');
          $wizardStepsItems.eq($(prevID).index() - 1).addClass('u-checked');
          $wizardStepsItems.eq($(prevID).index()).addClass('active');

          $wizardContentItems.hide().removeClass('active');
          $(prevID).fadeIn(400).addClass('active');
        });

        //Actions
        collection = collection.add($this);
      });
    }

  };

})(jQuery);
